CXX := clang++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -pedantic -fopenmp
LDFLAGS := -fopenmp
TARGET := dem_solver
SRC := src/main.cpp

.PHONY: all clean run-free-fall run-constant-velocity run-bounce run-verify run-experiment run-scaling run-neighbor run-science full

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

run-free-fall: $(TARGET)
	.\$(TARGET) free_fall results

run-constant-velocity: $(TARGET)
	.\$(TARGET) constant_velocity results

run-bounce: $(TARGET)
	.\$(TARGET) bounce results

run-verify: $(TARGET)
	.\$(TARGET) verification results

run-experiment: $(TARGET)
	.\$(TARGET) experiment results

run-scaling: $(TARGET)
	.\$(TARGET) scaling results

run-neighbor: $(TARGET)
	.\$(TARGET) neighbor_bonus results

run-science: $(TARGET)
	.\$(TARGET) science_bonus results

full:
	powershell -ExecutionPolicy Bypass -File .\run_all.ps1

clean:
	if exist $(TARGET).exe del /Q $(TARGET).exe
