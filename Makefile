CXX := g++
CXXFLAGS := -O3 -march=native -std=c++20 -pthread
BIN := bin

TARGETS := $(BIN)/campaign_2d $(BIN)/rpo_3d_fixedvz $(BIN)/trace_3d_family \
           $(BIN)/snap_3d_inertial $(BIN)/vertical_bifurcation $(BIN)/continuation \
           $(BIN)/render_gif

all: $(TARGETS)

$(BIN):
	mkdir -p $(BIN)

$(BIN)/%: src/%.cpp | $(BIN)
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -rf $(BIN)

.PHONY: all clean
