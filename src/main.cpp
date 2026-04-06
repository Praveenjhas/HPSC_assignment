#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

// Small 3D vector helper used for positions, velocities, and forces.
struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    //various operators defined for addition ,substraction,multiply of a vector
    Vec3 &operator+=(const Vec3 &rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    Vec3 &operator-=(const Vec3 &rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    Vec3 &operator*=(double s) {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }
};

Vec3 operator+(Vec3 lhs, const Vec3 &rhs) {
    lhs += rhs;
    return lhs;
}

Vec3 operator-(Vec3 lhs, const Vec3 &rhs) {
    lhs -= rhs;
    return lhs;
}

Vec3 operator*(Vec3 lhs, double s) {
    lhs *= s;
    return lhs;
}

Vec3 operator*(double s, Vec3 rhs) {
    rhs *= s;
    return rhs;
}

Vec3 operator/(Vec3 lhs, double s) {
    lhs.x /= s;
    lhs.y /= s;
    lhs.z /= s;
    return lhs;
}

double dot(const Vec3 &a, const Vec3 &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double norm(const Vec3 &v) {
    return std::sqrt(dot(v, v));
}

// Stores the state of one spherical particle.
struct Particle {
    Vec3 position;
    Vec3 velocity;
    Vec3 force;
    double radius = 0.01;//default radius taken 
    double mass = 1.0;//default mass taken
};

struct Domain {
    double lx = 1.0;
    double ly = 1.0;
    double lz = 1.0;
};

// Groups together the physical and numerical parameters for one simulation run.
struct SimulationParams {
    double dt = 1.0e-4; //time step for the simulation
    double total_time = 1.0;//total time of simulation so steps=total_time/dt
    double kn = 2.0e4; //spring stiffness
    double gamma_n = 20.0; //damping coefficient
    Vec3 gravity{0.0, 0.0, -9.81};  //gravity vectors
    int output_every = 20; //frequency of saving particle data i.e saving particle data every 20 steps
    int progress_every = 0;//controls how often progress is printed
    bool use_neighbor_search = false;//switch for spacial optimzation,check only nearby particles for collision(grid-based approach)
    bool use_openmp = false;//enable for the parallel computation
    bool clamp_contact_force = true;  //prevents negative force  if(f<0) f=0;
};

// Profiling buckets used for the report section on runtime distribution.
// Instead of an external profiling tool, the code uses std::chrono timers and
// accumulates how much wall-clock time is spent in each major stage.
struct ProfileData {
    double zero_force = 0.0;
    double gravity = 0.0;
    double particle_contacts = 0.0;//time spent during particle collisions
    double wall_contacts = 0.0;//time spent during 
    double integration = 0.0;//time spent updating positions and velocity
    double diagnostics = 0.0;//time spent  computing diagnostics quantity

    double total_contact_time() const {
        return particle_contacts + wall_contacts;  //total time spent in all contact computations
    }
};

// state summary of the system at any time t that are written to diagonistics.csv.
struct Diagnostics {
    double time = 0.0; //current simulation time=dt*step;
    double kinetic_energy = 0.0;//total energy in the system
    double max_speed = 0.0;//maximum velocity magnitude among all particles
    double center_of_mass_x = 0.0;  //com of the entire system coordinate, weighted by mass
    double center_of_mass_y = 0.0;
    double center_of_mass_z = 0.0;
    std::size_t active_contacts = 0;//Number of actual collision happening 
    std::size_t candidate_pairs = 0;//Number of pairs checked for the collision
    double max_height = 0.0;  //maximum z-coordinate among all particles
};

// Stores contact counter at one timestep such as active contacts of collision and candidates checked for the collision
struct StepStats {
    std::size_t active_contacts = 0;
    std::size_t candidate_pairs = 0;
};


// Final information returned after one full simulation run finishes.
// wall_seconds is the total runtime of the whole simulation, while profile
// stores the per-stage timing breakdown used in profiling and scaling studies.
struct RunSummary {
    double wall_seconds = 0.0;  //actual real cpu time  to run full one simulation measured using std::chrono
    double simulated_seconds = 0.0;//one simulation expected ideal time=steps*dt
    Diagnostics last;//the final state of the system at the end of the simulation
    ProfileData profile;//stores performance breakdown of the system in the simulation done
};


// Small std::chrono-based timer used to accumulate wall-clock time into one
// profiling bucket automatically for the lifetime of a code block.
class ScopedTimer {
public:
    explicit ScopedTimer(double &bucket) : bucket_(bucket), start_(Clock::now()) {}

    ~ScopedTimer() {
        const auto end = Clock::now();
        bucket_ += std::chrono::duration<double>(end - start_).count();
    }

private:
    using Clock = std::chrono::steady_clock;
    double &bucket_;
    Clock::time_point start_;
};

// Stores the input parameters of simulation ,which experiment is running free_fall,bounce etc and corresponding dt,total_time parameters
struct ConfigSnapshot {
    std::string mode = "experiment";
    std::map<std::string, std::string> values;
};

//converts double->string with high precision which is required for hpsc
std::string to_string_precise(double value) {
    std::ostringstream out;
    out << std::setprecision(16) << value;
    return out.str();
}

// Creates the output folder if it does not already exist.
void ensure_directory(const fs::path &path) {
    fs::create_directories(path);
}

// Writes one comma-separated row in csv file  used in diagnostics.csv ,particle outputs,experiments results
void write_csv_line(std::ofstream &out, const std::vector<std::string> &fields) {
    for (std::size_t i = 0; i < fields.size(); ++i) {
        if (i > 0) {
            out << ',';
        }
        out << fields[i];
    }
    out << '\n';
}

// Computes the normal contact force between a particle pair using the spring-dashpot model.
Vec3 compute_contact_force(const Particle &a, const Particle &b, double kn, double gamma_n, bool clamp_force) {
    const Vec3 rij = b.position - a.position;
    const double dist = norm(rij);
    if (dist <= 1.0e-12) {
        return {};
    }
    const double overlap = a.radius + b.radius - dist;
    if (overlap <= 0.0) {
        return {};
    }

    const Vec3 nij = rij / dist;
    const Vec3 rel_vel = b.velocity - a.velocity;
    const double vn = dot(rel_vel, nij);
    double fn = kn * overlap - gamma_n * vn;
    if (clamp_force && fn < 0.0) {
        fn = 0.0;
    }
    return fn * nij;
}

// Resets all particle forces before the next timestep starts.
void zero_forces(std::vector<Particle> &particles) {
    for (auto &p : particles) {
        p.force = {};
    }
}

// Adds the body force due to gravity.
void add_gravity(std::vector<Particle> &particles, const Vec3 &g) {
    for (auto &p : particles) {
        p.force += p.mass * g;
    }
}

// Computes one wall-contact contribution using the same normal model as particle contact.
Vec3 wall_force_from_overlap(double overlap, const Vec3 &normal, const Vec3 &velocity, double kn, double gamma_n, bool clamp_force) {
    if (overlap <= 0.0) {
        return {};
    }

    const double vn = dot(velocity, normal);
    double fn = kn * overlap - gamma_n * vn;
    if (clamp_force && fn < 0.0) {
        fn = 0.0;
    }
    return fn * normal;
}

// Checks all six box walls and adds wall-contact forces when overlap exists.
std::size_t compute_wall_contacts(std::vector<Particle> &particles, const Domain &box, const SimulationParams &params) {
    std::size_t contacts = 0;
    for (auto &p : particles) {
        if (const double overlap = p.radius - p.position.x; overlap > 0.0) {
            p.force += wall_force_from_overlap(overlap, {1.0, 0.0, 0.0}, p.velocity, params.kn, params.gamma_n, params.clamp_contact_force);
            ++contacts;
        }
        if (const double overlap = p.position.x + p.radius - box.lx; overlap > 0.0) {
            p.force += wall_force_from_overlap(overlap, {-1.0, 0.0, 0.0}, p.velocity, params.kn, params.gamma_n, params.clamp_contact_force);
            ++contacts;
        }
        if (const double overlap = p.radius - p.position.y; overlap > 0.0) {
            p.force += wall_force_from_overlap(overlap, {0.0, 1.0, 0.0}, p.velocity, params.kn, params.gamma_n, params.clamp_contact_force);
            ++contacts;
        }
        if (const double overlap = p.position.y + p.radius - box.ly; overlap > 0.0) {
            p.force += wall_force_from_overlap(overlap, {0.0, -1.0, 0.0}, p.velocity, params.kn, params.gamma_n, params.clamp_contact_force);
            ++contacts;
        }
        if (const double overlap = p.radius - p.position.z; overlap > 0.0) {
            p.force += wall_force_from_overlap(overlap, {0.0, 0.0, 1.0}, p.velocity, params.kn, params.gamma_n, params.clamp_contact_force);
            ++contacts;
        }
        if (const double overlap = p.position.z + p.radius - box.lz; overlap > 0.0) {
            p.force += wall_force_from_overlap(overlap, {0.0, 0.0, -1.0}, p.velocity, params.kn, params.gamma_n, params.clamp_contact_force);
            ++contacts;
        }
    }
    return contacts;
}

 
//we need optimization so that a particle only collides with nearby particles only ,not all
//so we divide space into grid cells and only check same cell and neighbouring cells for collision


// Integer grid cell locations used by the neighbor grid.
struct CellIndex {
    int ix = 0;   //cell x index
    int iy = 0;   //cell y index
    int iz = 0;     //cell y index
};

// Uniform Cartesian grid used to reduce the cost of contact collision particles searching.
class NeighborGrid {
public:
    // Builds a grid that covers the full box with a chosen cell size.
    NeighborGrid(const Domain &box, double cell_size)
        : box_(box), cell_size_(cell_size) {
        nx_ = std::max(1, static_cast<int>(std::ceil(box_.lx / cell_size_)));
        ny_ = std::max(1, static_cast<int>(std::ceil(box_.ly / cell_size_)));
        nz_ = std::max(1, static_cast<int>(std::ceil(box_.lz / cell_size_)));
        cells_.assign(static_cast<std::size_t>(nx_ * ny_ * nz_), {});
    }

    // Reassigns every particle to the correct grid cell for the current timestep.
    void rebuild(const std::vector<Particle> &particles) {
        for (auto &cell : cells_) {
            cell.clear();
        }
        for (std::size_t i = 0; i < particles.size(); ++i) {
            cells_[linear_index(index_for(particles[i].position))].push_back(static_cast<int>(i));
        }
    }

    // Loops over only the relevant neighboring cells and calls the callback for each candidate pair.
    template <typename Callback>
    std::size_t for_candidate_pairs(const std::vector<Particle> &particles, Callback callback) const {
        (void)particles;
        std::size_t candidate_pairs = 0;
        // Only same-cell and neighboring-cell pairs are checked, which cuts down the all-pairs cost.
        for (int iz = 0; iz < nz_; ++iz) {
            for (int iy = 0; iy < ny_; ++iy) {
                for (int ix = 0; ix < nx_; ++ix) {
                    const CellIndex cell{ix, iy, iz};
                    const auto &current = cells_[linear_index(cell)];
                    if (current.empty()) {
                        continue;
                    }

                    for (int dz = -1; dz <= 1; ++dz) {
                        for (int dy = -1; dy <= 1; ++dy) {
                            for (int dx = -1; dx <= 1; ++dx) {
                                const CellIndex other{ix + dx, iy + dy, iz + dz};
                                //skip invalid cells
                                if (!is_valid(other)) {
                                    continue;
                                }
                                //avoid double counting and ensures each pari processed once
                                if (lexicographically_before(other, cell)) {
                                    continue;
                                }

                                const auto &neighbor = cells_[linear_index(other)];
                                if (neighbor.empty()) {
                                    continue;
                                }

                                for (std::size_t a = 0; a < current.size(); ++a) {
                                    const std::size_t b_start = same_cell(cell, other) ? a + 1 : 0;
                                    for (std::size_t b = b_start; b < neighbor.size(); ++b) {
                                        ++candidate_pairs; //inncreasing the candidate pair for collisions count
                                        callback(current[a], neighbor[b]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return candidate_pairs;
    }

private:
    Domain box_;
    double cell_size_ = 0.1;
    int nx_ = 1;
    int ny_ = 1;
    int nz_ = 1;
    std::vector<std::vector<int>> cells_;
    //convert particle position->grid cell index
    CellIndex index_for(const Vec3 &pos) const {
     //ensure index stays inside grid
        auto clamp_index = [](int id, int upper) {
            return std::max(0, std::min(id, upper - 1));
        };
        return {
            clamp_index(static_cast<int>(pos.x / cell_size_), nx_),
            clamp_index(static_cast<int>(pos.y / cell_size_), ny_),
            clamp_index(static_cast<int>(pos.z / cell_size_), nz_)
        };
    }
   //convert 3d to 1d single index
    std::size_t linear_index(const CellIndex &idx) const {
        return static_cast<std::size_t>(idx.ix + nx_ * (idx.iy + ny_ * idx.iz));
    }
    //avoid invalid cells
    bool is_valid(const CellIndex &idx) const {
        return idx.ix >= 0 && idx.ix < nx_ && idx.iy >= 0 && idx.iy < ny_ && idx.iz >= 0 && idx.iz < nz_;
    }
      //avoid same cell 
    static bool same_cell(const CellIndex &a, const CellIndex &b) {
        return a.ix == b.ix && a.iy == b.iy && a.iz == b.iz;
    }
    //avoid duplicates
    static bool lexicographically_before(const CellIndex &a, const CellIndex &b) {
        return std::tie(a.iz, a.iy, a.ix) < std::tie(b.iz, b.iy, b.ix);
    }
};
 
//for all relevent particle pairs compute contact force,apply newton's 3rd law and count statistics
StepStats compute_particle_contacts_serial(std::vector<Particle> &particles, const SimulationParams &params, const std::optional<NeighborGrid> &grid) {
    // Serial contact routine used for the baseline solver and for the neighbor-grid version.
    StepStats stats;
     //if optimization is enabled used the grid  based searching else the brute force searching method
    if (grid.has_value()) {
        stats.candidate_pairs = grid->for_candidate_pairs(particles, [&](int ia, int ib) {
            Vec3 force = compute_contact_force(particles[ia], particles[ib], params.kn, params.gamma_n, params.clamp_contact_force);
            //checking if collision occurs by detecting overlap and applying force
            if (norm(force) > 0.0) {
                particles[ia].force -= force;
                particles[ib].force += force;
                ++stats.active_contacts;
            }
        });
        return stats;
    }
     
    const std::size_t n = particles.size();
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            ++stats.candidate_pairs;
            Vec3 force = compute_contact_force(particles[i], particles[j], params.kn, params.gamma_n, params.clamp_contact_force);
            if (norm(force) > 0.0) {
                particles[i].force -= force;
                particles[j].force += force;
                ++stats.active_contacts;
            }
        }
    }
    return stats;
}

 
// OpenMp Parallel brute-force collision detection using thread-local force buffers to avoid race conditions
StepStats compute_particle_contacts_parallel(std::vector<Particle> &particles, const SimulationParams &params) {
    // OpenMP version of the brute-force pair loop.
    StepStats stats;
#ifdef _OPENMP
    const std::size_t n = particles.size();
    const int threads = omp_get_max_threads();
     //Race conditions avoidence hence to define private varaibles
    // Each thread stores force contributions locally first, which avoids races on particle.force.
    std::vector<std::vector<Vec3>> force_buffers(static_cast<std::size_t>(threads), std::vector<Vec3>(n));
     //each thread tracks its own contact count ,its own pair count to avoid the contention
    std::vector<std::size_t> contacts(static_cast<std::size_t>(threads), 0);
    std::vector<std::size_t> candidates(static_cast<std::size_t>(threads), 0);

#pragma omp parallel
    { 
       
        const int tid = omp_get_thread_num();
        auto &local_forces = force_buffers[static_cast<std::size_t>(tid)];
        std::size_t local_contacts = 0;
        std::size_t local_candidates = 0;
//openmp splits loop among threads
#pragma omp for schedule(static)
        for (int i = 0; i < static_cast<int>(n); ++i) {
            for (std::size_t j = static_cast<std::size_t>(i) + 1; j < n; ++j) {
                ++local_candidates;
                Vec3 force = compute_contact_force(particles[static_cast<std::size_t>(i)], particles[j], params.kn, params.gamma_n, params.clamp_contact_force);
                if (norm(force) > 0.0) {
                    local_forces[static_cast<std::size_t>(i)] -= force;
                    local_forces[j] += force;
                    ++local_contacts;
                }
            }
        }

        contacts[static_cast<std::size_t>(tid)] = local_contacts;
        candidates[static_cast<std::size_t>(tid)] = local_candidates;
    }
     
   //combine contributions from all  threads
    for (const auto &buffer : force_buffers) {
        for (std::size_t i = 0; i < n; ++i) {
            particles[i].force += buffer[i];
        }
    }

    stats.active_contacts = std::accumulate(contacts.begin(), contacts.end(), std::size_t{0});
    stats.candidate_pairs = std::accumulate(candidates.begin(), candidates.end(), std::size_t{0});
#else
    stats = compute_particle_contacts_serial(particles, params, std::nullopt);
#endif
    return stats;
}

//update particle velocity and position using semi-implicit euler update
void integrate(std::vector<Particle> &particles, double dt) {
    // Semi-implicit Euler update: first velocity, then position.
    for (auto &p : particles) {
        p.velocity += (p.force / p.mass) * dt;
        p.position += p.velocity * dt;
    }
}
//computes summary of system at time t and dianostics data are calculated accordingly
Diagnostics compute_diagnostics(const std::vector<Particle> &particles, double time, const StepStats &step_stats) {
    // Collects the main scalar quantities used in plots, progress printing, and summaries.
    Diagnostics d;
    d.time = time;
    d.active_contacts = step_stats.active_contacts;
    d.candidate_pairs = step_stats.candidate_pairs;

    double mass_sum = 0.0;
    // Diagnostics are aggregated here so the same data can be used both for CSV output and live progress printing.
    for (const auto &p : particles) {
        const double speed = norm(p.velocity);
        d.kinetic_energy += 0.5 * p.mass * speed * speed;
        d.max_speed = std::max(d.max_speed, speed);
        d.center_of_mass_x += p.mass * p.position.x;
        d.center_of_mass_y += p.mass * p.position.y;
        d.center_of_mass_z += p.mass * p.position.z;
        d.max_height = std::max(d.max_height, p.position.z);
        mass_sum += p.mass;
    }
    if (mass_sum > 0.0) {
        d.center_of_mass_x /= mass_sum;
        d.center_of_mass_y /= mass_sum;
        d.center_of_mass_z /= mass_sum;
    }
    return d;
}



void write_particles_csv(const fs::path &path, const std::vector<Particle> &particles) {
    // Writes one particle snapshot so configurations can be plotted later.
    std::ofstream out(path);
    write_csv_line(out, {"id", "x", "y", "z", "vx", "vy", "vz", "radius", "mass"});
    for (std::size_t i = 0; i < particles.size(); ++i) {
        const auto &p = particles[i];
        write_csv_line(out, {
            std::to_string(i),
            to_string_precise(p.position.x),
            to_string_precise(p.position.y),
            to_string_precise(p.position.z),
            to_string_precise(p.velocity.x),
            to_string_precise(p.velocity.y),
            to_string_precise(p.velocity.z),
            to_string_precise(p.radius),
            to_string_precise(p.mass)
        });
    }
}

void write_config_file(const fs::path &path, const ConfigSnapshot &config) {
    // Stores the run settings in plain text for reproducibility.
    std::ofstream out(path);
    out << "mode=" << config.mode << '\n';
    for (const auto &[key, value] : config.values) {
        out << key << '=' << value << '\n';
    }
}

void write_diagnostics_header(std::ofstream &out) {
    // Header for the main timestep diagnostics file.
    write_csv_line(out, {
        "time", "kinetic_energy", "max_speed", "center_of_mass_x", "center_of_mass_y", "center_of_mass_z",
        "active_contacts", "candidate_pairs", "max_height"
    });
}

void append_diagnostics(std::ofstream &out, const Diagnostics &d) {
    // Appends one timestep worth of diagnostics to diagnostics.csv.
    write_csv_line(out, {
        to_string_precise(d.time),
        to_string_precise(d.kinetic_energy),
        to_string_precise(d.max_speed),
        to_string_precise(d.center_of_mass_x),
        to_string_precise(d.center_of_mass_y),
        to_string_precise(d.center_of_mass_z),
        std::to_string(d.active_contacts),
        std::to_string(d.candidate_pairs),
        to_string_precise(d.max_height)
    });
}

void print_run_summary(const RunSummary &summary, const std::string &label) {
    // Compact summary printed after one run finishes.
    std::cout << "\n[" << label << "]\n";
    std::cout << "simulated time     : " << summary.simulated_seconds << " s\n";
    std::cout << "wall-clock runtime : " << summary.wall_seconds << " s\n";
    std::cout << "final KE           : " << summary.last.kinetic_energy << '\n';
    std::cout << "max speed          : " << summary.last.max_speed << '\n';
    std::cout << "contacts           : " << summary.last.active_contacts << '\n';
    std::cout << "candidate pairs    : " << summary.last.candidate_pairs << '\n';
    std::cout << "profile contacts   : " << summary.profile.total_contact_time() << " s\n";
}

void print_run_banner(const std::string &label) {
    // Short label shown before the progress lines of each case.
    std::cout << "\n>>> Running " << label << '\n';
}

void print_progress_line(const Diagnostics &d, int step, int total_steps) {
    // Progress output shown during longer runs so the simulation state is visible in the terminal.
    std::cout << "step " << step << "/" << total_steps
              << " | t = " << std::fixed << std::setprecision(4) << d.time << " s"
              << " | KE = " << std::scientific << std::setprecision(4) << d.kinetic_energy
              << " | COM = (" << std::fixed << std::setprecision(4)
              << d.center_of_mass_x << ", " << d.center_of_mass_y << ", " << d.center_of_mass_z << ")"
              << " | contacts = " << d.active_contacts
              << " | candidate_pairs = " << d.candidate_pairs
              << " | max_speed = " << std::fixed << std::setprecision(4) << d.max_speed
              << '\n';
}

// Runs one complete simulation from start to finish.
// This is also where the assignment's profiling information is collected:
// total runtime, runtime distribution across major functions, and final diagnostics.
RunSummary run_simulation(std::vector<Particle> particles,
                          const Domain &box,
                          const SimulationParams &params,
                          const fs::path &output_dir,
                          const ConfigSnapshot &config,
                          bool write_particle_snapshots) {
    // Shared driver used by all assignment cases: setup, timestep loop, diagnostics, and summaries.
    ensure_directory(output_dir);
    write_config_file(output_dir / "run_config.txt", config);

    std::ofstream diag_out(output_dir / "diagnostics.csv");
    write_diagnostics_header(diag_out);

    std::optional<NeighborGrid> grid;
    if (params.use_neighbor_search) {
        double max_radius = 0.0;
        for (const auto &p : particles) {
            max_radius = std::max(max_radius, p.radius);
        }
        grid.emplace(box, std::max(2.2 * max_radius, 1.0e-3));
    }
    
    // Profiling data for one run.
    // The compiler optimization part of the assignment is handled at build time
    // with flags such as -O3 in the Makefile / build commands, not here at runtime.
    ProfileData profile;
    const int steps = static_cast<int>(std::ceil(params.total_time / params.dt));
    const int progress_stride = (params.progress_every > 0) ? params.progress_every : std::max(1, steps / 10);
    // Total runtime timer for the full simulation.
    auto wall_start = std::chrono::steady_clock::now();
    Diagnostics last;

    // Write the initial state before any integration so timestamps stay aligned with the physical state.
    last = compute_diagnostics(particles, 0.0, {});
    append_diagnostics(diag_out, last);
    print_progress_line(last, 0, steps);
    if (write_particle_snapshots) {
        std::ostringstream name;
        name << "particles_step_" << std::setw(6) << std::setfill('0') << 0 << ".csv";
        write_particles_csv(output_dir / name.str(), particles);
    }

    // This is the main timestep loop used by all experiments in the assignment.
    for (int step = 1; step <= steps; ++step) {
        const double time = step * params.dt;
        {
            // Stage timer: reset old forces.
            ScopedTimer timer(profile.zero_force);  //clears old forces
            zero_forces(particles);
        }
        {
            // Stage timer: apply gravity force.
            ScopedTimer timer(profile.gravity);  //add gravity
            add_gravity(particles, params.gravity);
        }

        if (params.use_neighbor_search && grid.has_value()) {   //reassign particles to cells
            grid->rebuild(particles);
        }

        StepStats stats;
        {
            // Stage timer: particle-particle contact work.
            // This bucket is usually the dominant cost and is the main profiling hotspot.
            ScopedTimer timer(profile.particle_contacts);
            // The neighbor-search path and the OpenMP path are kept separate here to keep each version easier to reason about.
            if (params.use_openmp && !params.use_neighbor_search) {
                stats = compute_particle_contacts_parallel(particles, params);
            } else {
                stats = compute_particle_contacts_serial(particles, params, grid);
            }
        }
        // adds wall contacts
        {
            // Stage timer: wall-contact handling.
            ScopedTimer timer(profile.wall_contacts);
            stats.active_contacts += compute_wall_contacts(particles, box, params);
        }
            //motion update
        {
            // Stage timer: time integration update.
            ScopedTimer timer(profile.integration);
            integrate(particles, params.dt);
        }
           
            //computing diagnositcs such as K.E,COM and contacts etc.
        {
            // Stage timer: diagnostics/output quantities used for plots and summaries.
            ScopedTimer timer(profile.diagnostics);
            last = compute_diagnostics(particles, time, stats);
            append_diagnostics(diag_out, last);
        }
       //printing to terminal 
        if (step == steps || step % progress_stride == 0) {
            print_progress_line(last, step, steps);
        }
     //saving particle data
        if (write_particle_snapshots && (step % params.output_every == 0 || step == steps)) {
            std::ostringstream name;
            name << "particles_step_" << std::setw(6) << std::setfill('0') << step << ".csv";
            write_particles_csv(output_dir / name.str(), particles);
        }
    }

    auto wall_end = std::chrono::steady_clock::now();
    // Total runtime and runtime distribution are packed into RunSummary here.
    // Later, summary.txt and CSV outputs use these values for profiling plots,
    // speedup/efficiency tables, and bottleneck discussion in the report.
    RunSummary summary;
    summary.wall_seconds = std::chrono::duration<double>(wall_end - wall_start).count();
    summary.simulated_seconds = steps * params.dt;
    summary.last = last;
    summary.profile = profile;
    return summary;
}

double sphere_mass(double radius, double density) {
    // Mass of one spherical particle assuming constant density.
    constexpr double pi = 3.14159265358979323846;
    return density * (4.0 / 3.0) * pi * radius * radius * radius;
}

//making a particle by its parameters
std::vector<Particle> make_single_particle(const Vec3 &pos, const Vec3 &vel, double radius, double density) {
    // Helper for one-particle test cases.
    Particle p;
    p.position = pos;
    p.velocity = vel;
    p.radius = radius;
    p.mass = sphere_mass(radius, density);
    return {p};
}

//creating a random cloud of particle inside a 3D box 
std::vector<Particle> make_particle_cloud(int count, const Domain &box, double radius, double density, unsigned seed) {
    // Creates a random particle cloud for multi-particle experiments.
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> xdist(radius, box.lx - radius);
    std::uniform_real_distribution<double> ydist(radius, box.ly - radius);
    std::uniform_real_distribution<double> zdist(box.lz * 0.45, box.lz - radius);
    std::uniform_real_distribution<double> jitter(-0.05, 0.05);

    std::vector<Particle> particles;
    particles.reserve(static_cast<std::size_t>(count));

    for (int i = 0; i < count; ++i) {
        Particle p;
        // A small random initial velocity helps break perfectly symmetric starts in cloud-settling cases.
        p.position = {xdist(rng), ydist(rng), zdist(rng)};
        p.velocity = {jitter(rng), jitter(rng), jitter(rng)};
        p.radius = radius;
        p.mass = sphere_mass(radius, density);
        particles.push_back(p);
    }
    return particles;
}

std::vector<Particle> make_particle_cloud_zero_velocity(int count, const Domain &box, double radius, double density, unsigned seed) {
    // Cloud used for settling studies where the assignment asks for zero initial velocity.
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> xdist(radius, box.lx - radius);
    std::uniform_real_distribution<double> ydist(radius, box.ly - radius);
    std::uniform_real_distribution<double> zdist(box.lz * 0.45, box.lz - radius);

    std::vector<Particle> particles;
    particles.reserve(static_cast<std::size_t>(count));

    for (int i = 0; i < count; ++i) {
        Particle p;
        p.position = {xdist(rng), ydist(rng), zdist(rng)};
        p.velocity = {0.0, 0.0, 0.0};
        p.radius = radius;
        p.mass = sphere_mass(radius, density);
        particles.push_back(p);
    }
    return particles;
}

//readble summary of simulation results
void write_text_summary(const fs::path &path, const std::vector<std::pair<std::string, std::string>> &lines) {
    // Writes a short plain-text summary file for quick inspection.
    std::ofstream out(path);
    for (const auto &[key, value] : lines) {
        out << std::left << std::setw(28) << key << ": " << value << '\n';
    }
}

// Verification case 1: compare one-particle free fall against the analytical solution before floor impact.
void run_free_fall_case(const fs::path &root) {
    print_run_banner("free_fall");
    const Domain box{1.0, 1.0, 3.0};
    SimulationParams params;
    params.dt = 1.0e-4;
    // End this case before first wall contact so the analytical free-fall comparison stays valid.
    params.total_time = 0.6;
    params.gravity = {0.0, 0.0, -9.81};
    params.kn = 2.0e5;
    params.gamma_n = 25.0;
    params.output_every = 400;

    const double z0 = 2.0;
    const auto particles = make_single_particle({0.5, 0.5, z0}, {0.0, 0.0, 0.0}, 0.03, 2500.0);
    const fs::path out_dir = root / "free_fall";

    ConfigSnapshot cfg;
    cfg.mode = "free_fall";
    cfg.values = {
        {"dt", to_string_precise(params.dt)},
        {"total_time", to_string_precise(params.total_time)},
        {"z0", to_string_precise(z0)}
    };

    const RunSummary summary = run_simulation(particles, box, params, out_dir, cfg, true);
    const double final_z_exact = z0 + 0.5 * params.gravity.z * params.total_time * params.total_time;
    const double final_v_exact = params.gravity.z * params.total_time;

    write_text_summary(out_dir / "summary.txt", {
        {"final_numeric_ke", to_string_precise(summary.last.kinetic_energy)},
        {"final_exact_z", to_string_precise(final_z_exact)},
        {"final_exact_vz", to_string_precise(final_v_exact)},
        {"wall_runtime", to_string_precise(summary.wall_seconds)}
    });
    print_run_summary(summary, "free_fall");
}

// Verification case 2: with gravity disabled, the particle should keep a constant velocity in all directions.
void run_constant_velocity_case(const fs::path &root) {
    print_run_banner("constant_velocity");
    const Domain box{2.0, 2.0, 2.0};
    SimulationParams params;
    params.dt = 2.0e-4;
    params.total_time = 0.5;
    params.gravity = {0.0, 0.0, 0.0};
    params.output_every = 500;

    const Vec3 start{0.5, 0.6, 0.7};
    const Vec3 velocity{0.8, -0.35, 0.5};
    const auto particles = make_single_particle(start, velocity, 0.03, 2500.0);
    const fs::path out_dir = root / "constant_velocity";

    ConfigSnapshot cfg;
    cfg.mode = "constant_velocity";
    cfg.values = {
        {"dt", to_string_precise(params.dt)},
        {"total_time", to_string_precise(params.total_time)}
    };

    const RunSummary summary = run_simulation(particles, box, params, out_dir, cfg, true);
    const Vec3 exact = start + velocity * params.total_time;

    write_text_summary(out_dir / "summary.txt", {
        {"exact_x", to_string_precise(exact.x)},
        {"exact_y", to_string_precise(exact.y)},
        {"exact_z", to_string_precise(exact.z)},
        {"runtime", to_string_precise(summary.wall_seconds)}
    });
    print_run_summary(summary, "constant_velocity");
}

// Verification case 3: drop one particle onto the floor to check bounded wall penetration and damping-dependent rebound.
void run_bounce_case(const fs::path &root, double gamma_n, const std::string &label) {
    print_run_banner(label);
    const Domain box{1.0, 1.0, 1.2};
    SimulationParams params;
    params.dt = 5.0e-5;
    params.total_time = 2.0;
    params.gravity = {0.0, 0.0, -9.81};
    params.kn = 4.0e5;
    params.gamma_n = gamma_n;
    params.output_every = 2000;

    const auto particles = make_single_particle({0.5, 0.5, 0.8}, {0.0, 0.0, 0.0}, 0.04, 2500.0);
    const fs::path out_dir = root / label;

    ConfigSnapshot cfg;
    cfg.mode = "bounce";
    cfg.values = {
        {"gamma_n", to_string_precise(gamma_n)},
        {"dt", to_string_precise(params.dt)}
    };

    const RunSummary summary = run_simulation(particles, box, params, out_dir, cfg, true);
    write_text_summary(out_dir / "summary.txt", {
        {"gamma_n", to_string_precise(gamma_n)},
        {"runtime", to_string_precise(summary.wall_seconds)},
        {"final_ke", to_string_precise(summary.last.kinetic_energy)}
    });
    print_run_summary(summary, label);
}


// Runs a performance benchmark by simulating a particle system with a specified number of particles.
// It allows comparison between brute-force and optimized (neighbor grid / OpenMP) approaches by
// measuring runtime, contact computation cost, and number of candidate pairs. The results are saved
// for analysis of scalability and efficiency.
void run_particle_count_experiment(const fs::path &root, int count, bool use_openmp, bool use_neighbor_search, int threads) {
    // Multi-particle benchmark used for profiling and runtime comparisons.
    const Domain box{1.0, 1.0, 1.0};
    SimulationParams params;
    params.dt = 1.0e-4;
    params.total_time = 0.03;
    params.gravity = {0.0, 0.0, -9.81};
    params.kn = 2.0e4;
    params.gamma_n = 30.0;
    params.output_every = 500;
    params.use_openmp = use_openmp;
    params.use_neighbor_search = use_neighbor_search;

#ifdef _OPENMP
    if (use_openmp && threads > 0) {
        omp_set_num_threads(threads);
    }
#else
    (void)threads;
#endif

    const auto particles = make_particle_cloud(count, box, 0.015, 2400.0, static_cast<unsigned>(1234 + count));
    std::ostringstream name;
    name << "n" << count << (use_openmp ? "_omp" : "_serial") << (use_neighbor_search ? "_grid" : "_pairs");
    if (use_openmp) {
        name << "_t" << threads;
    }
    print_run_banner(name.str());
    const fs::path out_dir = root / name.str();

    ConfigSnapshot cfg;
    cfg.mode = "particle_count_experiment";
    cfg.values = {
        {"count", std::to_string(count)},
        {"use_openmp", use_openmp ? "1" : "0"},
        {"use_neighbor_search", use_neighbor_search ? "1" : "0"},
        {"threads", std::to_string(threads)}
    };

    const RunSummary summary = run_simulation(particles, box, params, out_dir, cfg, false);
    write_text_summary(out_dir / "summary.txt", {
        {"count", std::to_string(count)},
        {"runtime", to_string_precise(summary.wall_seconds)},
        {"zero_force_time", to_string_precise(summary.profile.zero_force)},
        {"gravity_time", to_string_precise(summary.profile.gravity)},
        {"particle_contact_time", to_string_precise(summary.profile.particle_contacts)},
        {"wall_contact_time", to_string_precise(summary.profile.wall_contacts)},
        {"integration_time", to_string_precise(summary.profile.integration)},
        {"diagnostics_time", to_string_precise(summary.profile.diagnostics)},
        {"contact_time", to_string_precise(summary.profile.total_contact_time())},
        {"candidate_pairs", std::to_string(summary.last.candidate_pairs)}
    });
    print_run_summary(summary, name.str());
}

void run_custom_particle_case(const fs::path &root, int count, double dt, double total_time) {
    // One-off manual timing case for arbitrary particle count and simulation length.
    const Domain box{1.0, 1.0, 1.0};
    SimulationParams params;
    params.dt = dt;
    params.total_time = total_time;
    params.gravity = {0.0, 0.0, -9.81};
    params.kn = 2.0e4;
    params.gamma_n = 30.0;
    params.output_every = std::numeric_limits<int>::max();
    params.use_openmp = false;
    params.use_neighbor_search = false;

    const auto particles = make_particle_cloud(count, box, 0.015, 2400.0, static_cast<unsigned>(7000 + count));

    std::ostringstream name;
    name << "custom_n" << count << "_dt" << std::scientific << dt << "_T" << total_time;
    print_run_banner(name.str());

    ConfigSnapshot cfg;
    cfg.mode = "custom_particle_case";
    cfg.values = {
        {"count", std::to_string(count)},
        {"dt", to_string_precise(dt)},
        {"total_time", to_string_precise(total_time)}
    };

    const fs::path out_dir = root / name.str();
    const RunSummary summary = run_simulation(particles, box, params, out_dir, cfg, false);
    write_text_summary(out_dir / "summary.txt", {
        {"count", std::to_string(count)},
        {"dt", to_string_precise(dt)},
        {"total_time", to_string_precise(total_time)},
        {"steps", std::to_string(static_cast<int>(std::ceil(total_time / dt)))},
        {"runtime", to_string_precise(summary.wall_seconds)},
        {"particle_contact_time", to_string_precise(summary.profile.particle_contacts)},
        {"diagnostics_time", to_string_precise(summary.profile.diagnostics)},
        {"candidate_pairs", std::to_string(summary.last.candidate_pairs)}
    });
    print_run_summary(summary, name.str());
}

// Evaluates parallel scalability by running the simulation with different OpenMP thread counts.
// Computes runtime, speedup (relative to single-thread), and efficiency to analyze how well
// the algorithm utilizes multiple CPU cores.
void run_thread_scaling(const fs::path &root, int count, const std::vector<int> &threads) {
    // Repeats the same case with different OpenMP thread counts to measure speedup and efficiency.
    ensure_directory(root);
    std::ofstream table(root / "speedup.csv");
    write_csv_line(table, {"threads", "runtime", "speedup", "efficiency"});

    double baseline = 0.0;
    // Thread counts are benchmarked one after another so the measured runtimes are comparable.
    for (int t : threads) {
        const Domain box{1.0, 1.0, 1.0};
        SimulationParams params;
        params.dt = 1.0e-4;
        params.total_time = 0.02;
        params.gravity = {0.0, 0.0, -9.81};
        params.kn = 2.0e4;
        params.gamma_n = 30.0;
        params.use_openmp = (t > 1);

#ifdef _OPENMP
        omp_set_num_threads(t);
#endif

        const auto particles = make_particle_cloud(count, box, 0.015, 2400.0, 777u);
        ConfigSnapshot cfg;
        cfg.mode = "thread_scaling";
        cfg.values = {
            {"count", std::to_string(count)},
            {"threads", std::to_string(t)}
        };

        std::ostringstream run_name;
        run_name << "threads_" << t;
        print_run_banner(run_name.str());
        const RunSummary summary = run_simulation(particles, box, params, root / run_name.str(), cfg, false);

        if (t == 1) {
            baseline = summary.wall_seconds;
        }

        const double speedup = baseline > 0.0 ? baseline / summary.wall_seconds : 1.0;
        const double efficiency = speedup / static_cast<double>(t);
        write_csv_line(table, {
            std::to_string(t),
            to_string_precise(summary.wall_seconds),
            to_string_precise(speedup),
            to_string_precise(efficiency)
        });
    }
}

// Measures weak scaling by increasing the particle count with the thread count
// so the approximate amount of work per thread stays comparable.
void run_weak_scaling(const fs::path &root, int base_count_per_thread, const std::vector<int> &threads) {
    ensure_directory(root);
    std::ofstream table(root / "weak_scaling.csv");
    write_csv_line(table, {"threads", "count", "box_length", "runtime", "weak_efficiency", "overhead"});

    double baseline = 0.0;
    for (int t : threads) {
        const int count = base_count_per_thread * t;
        // Increase the box volume with particle count so the packing density stays roughly constant.
        const double side = std::cbrt(static_cast<double>(t));
        const Domain box{side, side, side};
        SimulationParams params;
        params.dt = 1.0e-4;
        params.total_time = 0.02;
        params.gravity = {0.0, 0.0, -9.81};
        params.kn = 2.0e4;
        params.gamma_n = 30.0;
        params.use_openmp = (t > 1);

#ifdef _OPENMP
        omp_set_num_threads(t);
#endif

        const auto particles = make_particle_cloud(count, box, 0.015, 2400.0, static_cast<unsigned>(900 + t));
        ConfigSnapshot cfg;
        cfg.mode = "weak_scaling";
        cfg.values = {
            {"count", std::to_string(count)},
            {"threads", std::to_string(t)},
            {"base_count_per_thread", std::to_string(base_count_per_thread)},
            {"box_length", to_string_precise(side)}
        };

        std::ostringstream run_name;
        run_name << "threads_" << t << "_n" << count;
        print_run_banner(run_name.str());
        const RunSummary summary = run_simulation(particles, box, params, root / run_name.str(), cfg, false);

        if (t == 1) {
            baseline = summary.wall_seconds;
        }

        const double weak_efficiency = baseline > 0.0 ? baseline / summary.wall_seconds : 1.0;
        const double overhead = summary.wall_seconds - baseline;
        write_csv_line(table, {
            std::to_string(t),
            std::to_string(count),
            to_string_precise(side),
            to_string_precise(summary.wall_seconds),
            to_string_precise(weak_efficiency),
            to_string_precise(overhead)
        });
    }
}

// Verifies that the OpenMP version preserves numerical correctness by comparing
// a serial and parallel run for the same physical setup.
void run_parallel_correctness_check(const fs::path &root, int count, int threads) {
    ensure_directory(root);
    const bool append = fs::exists(root / "parallel_correctness.csv");

    const Domain box{1.0, 1.0, 1.0};
    SimulationParams serial_params;
    serial_params.dt = 1.0e-4;
    serial_params.total_time = 0.02;
    serial_params.gravity = {0.0, 0.0, -9.81};
    serial_params.kn = 2.0e4;
    serial_params.gamma_n = 30.0;
    serial_params.use_openmp = false;
    serial_params.output_every = 20;

    SimulationParams parallel_params = serial_params;
    parallel_params.use_openmp = true;

    const auto serial_particles = make_particle_cloud(count, box, 0.015, 2400.0, 2026u);
    const auto parallel_particles = make_particle_cloud(count, box, 0.015, 2400.0, 2026u);

    ConfigSnapshot serial_cfg;
    serial_cfg.mode = "parallel_correctness_serial";
    serial_cfg.values = {{"count", std::to_string(count)}};

    ConfigSnapshot parallel_cfg;
    parallel_cfg.mode = "parallel_correctness_parallel";
    parallel_cfg.values = {{"count", std::to_string(count)}, {"threads", std::to_string(threads)}};

#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    print_run_banner("parallel_correctness_serial");
    const fs::path serial_out = root / ("serial_reference_n" + std::to_string(count));
    const RunSummary serial_summary = run_simulation(serial_particles, box, serial_params, serial_out, serial_cfg, true);

#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    print_run_banner("parallel_correctness_parallel");
    const fs::path parallel_out = root / ("parallel_reference_n" + std::to_string(count));
    const RunSummary parallel_summary = run_simulation(parallel_particles, box, parallel_params, parallel_out, parallel_cfg, true);

    const double ke_abs_diff = std::abs(serial_summary.last.kinetic_energy - parallel_summary.last.kinetic_energy);
    const double max_speed_abs_diff = std::abs(serial_summary.last.max_speed - parallel_summary.last.max_speed);
    const double com_x_abs_diff = std::abs(serial_summary.last.center_of_mass_x - parallel_summary.last.center_of_mass_x);
    const double com_y_abs_diff = std::abs(serial_summary.last.center_of_mass_y - parallel_summary.last.center_of_mass_y);
    const double com_z_abs_diff = std::abs(serial_summary.last.center_of_mass_z - parallel_summary.last.center_of_mass_z);

    std::ofstream table;
    if (append) {
        table.open(root / "parallel_correctness.csv", std::ios::app);
    } else {
        table.open(root / "parallel_correctness.csv");
        write_csv_line(table, {
            "count",
            "threads",
            "serial_ke",
            "parallel_ke",
            "ke_abs_diff",
            "serial_max_speed",
            "parallel_max_speed",
            "max_speed_abs_diff",
            "serial_com_z",
            "parallel_com_z",
            "com_z_abs_diff",
            "serial_contacts",
            "parallel_contacts"
        });
    }
    write_csv_line(table, {
        std::to_string(count),
        std::to_string(threads),
        to_string_precise(serial_summary.last.kinetic_energy),
        to_string_precise(parallel_summary.last.kinetic_energy),
        to_string_precise(ke_abs_diff),
        to_string_precise(serial_summary.last.max_speed),
        to_string_precise(parallel_summary.last.max_speed),
        to_string_precise(max_speed_abs_diff),
        to_string_precise(serial_summary.last.center_of_mass_z),
        to_string_precise(parallel_summary.last.center_of_mass_z),
        to_string_precise(com_z_abs_diff),
        std::to_string(serial_summary.last.active_contacts),
        std::to_string(parallel_summary.last.active_contacts)
    });

    write_text_summary(root / ("summary_n" + std::to_string(count) + ".txt"), {
        {"count", std::to_string(count)},
        {"threads", std::to_string(threads)},
        {"ke_abs_diff", to_string_precise(ke_abs_diff)},
        {"max_speed_abs_diff", to_string_precise(max_speed_abs_diff)},
        {"com_x_abs_diff", to_string_precise(com_x_abs_diff)},
        {"com_y_abs_diff", to_string_precise(com_y_abs_diff)},
        {"com_z_abs_diff", to_string_precise(com_z_abs_diff)}
    });
}

void run_timestep_sensitivity(const fs::path &root) {
    // Studies how the timestep changes the free-fall error.
    ensure_directory(root);
    std::ofstream out(root / "timestep_scan.csv");
    write_csv_line(out, {"dt", "final_z_error", "runtime"});

    const std::vector<double> timesteps{1.0e-3, 5.0e-4, 2.5e-4, 1.0e-4};
    const Domain box{1.0, 1.0, 3.0};
    const double z0 = 2.0;

    for (double dt : timesteps) {
        SimulationParams params;
        params.dt = dt;
        params.total_time = 0.6;
        params.gravity = {0.0, 0.0, -9.81};
        params.kn = 2.0e5;
        params.gamma_n = 25.0;

        const auto particles = make_single_particle({0.5, 0.5, z0}, {0.0, 0.0, 0.0}, 0.03, 2500.0);
        std::ostringstream tag;
        tag << "dt_" << std::scientific << dt;

        ConfigSnapshot cfg;
        cfg.mode = "timestep_sensitivity";
        cfg.values = {{"dt", to_string_precise(dt)}};

        const RunSummary summary = run_simulation(particles, box, params, root / tag.str(), cfg, false);
        const double exact_z = z0 + 0.5 * params.gravity.z * params.total_time * params.total_time;
        const double final_z_estimate = summary.last.center_of_mass_z;

        write_csv_line(out, {
            to_string_precise(dt),
            to_string_precise(std::abs(final_z_estimate - exact_z)),
            to_string_precise(summary.wall_seconds)
        });
    }
}

void run_neighbor_search_benchmark(const fs::path &root) {
    // Compares the brute-force and neighbor-grid contact search strategies.
    ensure_directory(root);
    std::ofstream table(root / "neighbor_search_comparison.csv");
    write_csv_line(table, {
        "count", "mode", "runtime", "contact_stage_time", "candidate_pairs", "active_contacts"
    });

    const std::vector<int> counts{200, 1000, 5000};
    for (int count : counts) {
        for (bool use_grid : {false, true}) {
            const Domain box{1.0, 1.0, 1.0};
            SimulationParams params;
            params.dt = 1.0e-4;
            params.total_time = 0.01;
            params.gravity = {0.0, 0.0, -9.81};
            params.kn = 2.0e4;
            params.gamma_n = 25.0;
            params.use_neighbor_search = use_grid;

            const auto particles = make_particle_cloud(count, box, 0.012, 2400.0, static_cast<unsigned>(900 + count));

            ConfigSnapshot cfg;
            cfg.mode = use_grid ? "neighbor_grid" : "all_pairs";
            cfg.values = {{"count", std::to_string(count)}};

            const std::string run_label = std::string(use_grid ? "grid_" : "pairs_") + std::to_string(count);
            print_run_banner(run_label);

            const RunSummary summary = run_simulation(
                particles,
                box,
                params,
                root / run_label,
                cfg,
                false
            );

            write_csv_line(table, {
                std::to_string(count),
                use_grid ? "grid" : "all_pairs",
                to_string_precise(summary.wall_seconds),
                to_string_precise(summary.profile.particle_contacts),
                std::to_string(summary.last.candidate_pairs),
                std::to_string(summary.last.active_contacts)
            });
        }
    }
}

void run_damping_investigation(const fs::path &root) {
    // Bonus study: varies damping and records how rebound and energy change.
    const std::vector<double> gammas{5.0, 20.0, 50.0, 100.0};
    ensure_directory(root);
    std::ofstream out(root / "damping_scan.csv");
    write_csv_line(out, {"gamma_n", "runtime", "final_ke"});

    for (double gamma : gammas) {
        const std::string folder = "bounce_gamma_" + std::to_string(static_cast<int>(gamma));
        run_bounce_case(root, gamma, folder);

        const fs::path bounce_summary = root / folder / "summary.txt";
        std::ifstream in(bounce_summary);
        std::string line;
        double runtime = 0.0;
        double final_ke = 0.0;
        while (std::getline(in, line)) {
            const auto colon = line.find(std::string(":"));
            if (colon == std::string::npos) {
                continue;
            }
            const std::string key = line.substr(0, colon);
            const std::string value = line.substr(colon + 1);
            if (key.find("runtime") != std::string::npos) {
                runtime = std::stod(value);
            }
            if (key.find("final_ke") != std::string::npos) {
                final_ke = std::stod(value);
            }
        }
        write_csv_line(out, {
            to_string_precise(gamma),
            to_string_precise(runtime),
            to_string_precise(final_ke)
        });
    }

    const Domain box{1.0, 1.0, 1.0};
    std::ofstream cloud_out(root / "cloud_settling_by_damping.csv");
    write_csv_line(cloud_out, {"gamma_n", "runtime", "final_ke", "com_z"});
    for (double gamma : gammas) {
        SimulationParams params;
        params.dt = 1.0e-4;
        params.total_time = 0.03;
        params.gravity = {0.0, 0.0, -9.81};
        params.kn = 2.0e4;
        params.gamma_n = gamma;
        params.use_neighbor_search = true;

        const auto particles = make_particle_cloud_zero_velocity(1000, box, 0.012, 2400.0, static_cast<unsigned>(1000 + gamma));
        ConfigSnapshot cfg;
        cfg.mode = "cloud_settling_damping";
        cfg.values = {{"gamma_n", to_string_precise(gamma)}};

        const RunSummary summary = run_simulation(
            particles,
            box,
            params,
            root / ("cloud_gamma_" + std::to_string(static_cast<int>(gamma))),
            cfg,
            false
        );

        write_csv_line(cloud_out, {
            to_string_precise(gamma),
            to_string_precise(summary.wall_seconds),
            to_string_precise(summary.last.kinetic_energy),
            to_string_precise(summary.last.center_of_mass_z)
        });
    }
}

void run_cloud_settling(const fs::path &root) {
    print_run_banner("cloud_settling");
    // A small particle cloud settling under gravity gives a more system-level test than the one-particle cases.
    const Domain box{1.0, 1.0, 1.0};
    SimulationParams params;
    params.dt = 1.0e-4;
    params.total_time = 0.05;
    params.gravity = {0.0, 0.0, -9.81};
    params.kn = 2.0e4;
    params.gamma_n = 30.0;
    params.output_every = 125;
    params.use_neighbor_search = true;

    const auto particles = make_particle_cloud_zero_velocity(1200, box, 0.012, 2400.0, 42u);
    ConfigSnapshot cfg;
    cfg.mode = "particle_cloud_settling";
    cfg.values = {{"count", "1200"}};
    const RunSummary summary = run_simulation(particles, box, params, root / "cloud_settling", cfg, true);
    print_run_summary(summary, "cloud_settling");
}

void print_usage() {
    // Command-line help shown when no valid mode is given.
    std::cout
        << "Usage: dem_solver <mode> [output_dir]\n\n"
        << "Modes:\n"
        << "  free_fall\n"
        << "  constant_velocity\n"
        << "  bounce\n"
        << "  verification\n"
        << "  experiment\n"
        << "  scaling\n"
        << "  neighbor_bonus\n"
        << "  science_bonus\n"
        << "  custom_case <count> <dt> <total_time>\n";
}

int main(int argc, char **argv) {
    // Entry point: chooses which assignment case or grouped workflow to run.
    try {
        if (argc < 2) {
            print_usage();
            return 1;
        }

        const std::string mode = argv[1];
        const fs::path root = (argc >= 3) ? fs::path(argv[2]) : fs::path("results");

        // Each mode corresponds to one assignment task or one grouped set of runs used in the report.
        if (mode == "free_fall") {
            run_free_fall_case(root);
        } else if (mode == "constant_velocity") {
            run_constant_velocity_case(root);
        } else if (mode == "bounce") {
            run_bounce_case(root, 40.0, "bounce_default");
        } else if (mode == "verification") {
            run_free_fall_case(root / "verification");
            run_constant_velocity_case(root / "verification");
            run_bounce_case(root / "verification", 20.0, "bounce_light_damping");
            run_bounce_case(root / "verification", 70.0, "bounce_heavy_damping");
            run_timestep_sensitivity(root / "verification" / "timestep_sensitivity");
        } else if (mode == "experiment") {
            for (int count : {200, 1000, 5000}) {
                run_particle_count_experiment(root / "experiments", count, false, false, 1);
            }
        } else if (mode == "scaling") {
            for (int count : {200, 1000, 5000}) {
                run_particle_count_experiment(root / "profiling", count, false, false, 1);
            }
            // Strong scaling keeps the problem size fixed and varies only the thread count.
            run_thread_scaling(root / "openmp_scaling_n1000", 1000, {1, 2, 4, 8});
            run_thread_scaling(root / "openmp_scaling_n5000", 5000, {1, 2, 4, 8});
            // Weak scaling increases the particle count with the thread count.
            run_weak_scaling(root / "weak_scaling", 625, {1, 2, 4, 8});
            run_parallel_correctness_check(root / "parallel_correctness", 1000, 8);
            run_parallel_correctness_check(root / "parallel_correctness", 5000, 8);
        } else if (mode == "neighbor_bonus") {
            run_neighbor_search_benchmark(root / "neighbor_bonus");
        } else if (mode == "science_bonus") {
            run_damping_investigation(root / "science_bonus");
            run_cloud_settling(root / "science_bonus");
        } else if (mode == "custom_case") {
            if (argc < 5) {
                print_usage();
                return 1;
            }
            const int count = std::stoi(argv[2]);
            const double dt = std::stod(argv[3]);
            const double total_time = std::stod(argv[4]);
            const fs::path custom_root = (argc >= 6) ? fs::path(argv[5]) : fs::path("results");
            run_custom_particle_case(custom_root, count, dt, total_time);
        } else {
            print_usage();
            return 1;
        }

        return 0;
    } catch (const std::exception &ex) {
        std::cerr << "error: " << ex.what() << '\n';
        return 2;
    }
}
