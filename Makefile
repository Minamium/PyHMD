CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
TEST_OUTPUT ?= runs/verification
NATIVE := src/pymhd/native
BUILD_DIR ?= $(CURDIR)/build
FP_MODE ?= default
.PHONY: all cpu cpu-mpi cuda cuda-single test test-python
.PHONY: test-convergence test-gpu benchmark
CONVERGENCE_OUTPUT ?= runs/convergence
GPU_OUTPUT ?= runs/phase0_gpu
FROZEN_SOURCE ?=
test-convergence: cpu
	PYTHONPATH=src python3 tests/convergence/run.py --binary build/pymhd_cpu --output "$(CONVERGENCE_OUTPUT)"
test-gpu:
	@test -n "$(FROZEN_SOURCE)" || { echo 'Set FROZEN_SOURCE to an independently archived source directory with ../freeze.json'; exit 2; }
	PYTHONPATH=src python3 tests/phase0/gpu_suite.py --output "$(GPU_OUTPUT)" --reference-source "$(FROZEN_SOURCE)"
benchmark: cpu
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/cache_benchmark.cpp -o build/cache_benchmark
	build/cache_benchmark 600 128 3
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/contact_mixing_diagnostic.cpp -o build/contact_mixing_diagnostic
	build/contact_mixing_diagnostic
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/gpu_rhs_equivalence.cpp -o build/cached_rhs_cpu
	build/cached_rhs_cpu 64 32
all: cpu
cpu cpu-mpi cuda cuda-single:
	$(MAKE) -C $(NATIVE) $@ BUILD_DIR="$(abspath $(BUILD_DIR))"
test-python: cpu
	PYTHONPATH=src python3 -m unittest discover -s tests/python -v
test: cpu
	python3 tests/verify.py --binary build/pymhd_cpu --output "$(TEST_OUTPUT)"
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/reconstruction_equivalence.cpp -o build/reconstruction_equivalence
	build/reconstruction_equivalence
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/ct_invariants.cpp -o build/ct_invariants
	build/ct_invariants
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/paper_initial.cpp -o build/test_paper_initial
	build/test_paper_initial
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/pressure_recovery.cpp -o build/pressure_recovery
	build/pressure_recovery
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/thermal_recovery.cpp -o build/thermal_recovery
	build/thermal_recovery
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/thermal_robustness.cpp -o build/thermal_robustness
	build/thermal_robustness
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/step_retry.cpp -o build/step_retry
	build/step_retry
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/contact_advection.cpp -o build/contact_advection
	build/contact_advection
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/strong_contact.cpp -o build/strong_contact
	build/strong_contact
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/adiabatic_flow.cpp -o build/adiabatic_flow
	build/adiabatic_flow
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/boundary_scan.cpp -o build/boundary_scan
	build/boundary_scan
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/eta_cache.cpp -o build/eta_cache
	build/eta_cache
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/rate_reduction.cpp -o build/rate_reduction
	build/rate_reduction
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/checkpoint.cpp -o build/checkpoint
	build/checkpoint
	python3 tests/checkpoint_cli.py --binary build/pymhd_cpu
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/metric_reconstruction.cpp -o build/metric_reconstruction
	build/metric_reconstruction
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/resistive_geometry.cpp -o build/resistive_geometry
	build/resistive_geometry
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/spherical_resistive_operator.cpp -o build/spherical_resistive_operator
	build/spherical_resistive_operator
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/rejected_trial_diagnostics.cpp -o build/rejected_trial_diagnostics
	build/rejected_trial_diagnostics
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/eta_model_override.cpp -o build/eta_model_override
	build/eta_model_override
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/conservation_budget.cpp -o build/conservation_budget
	build/conservation_budget
	python3 tests/conservation_budget_cli.py --binary build/pymhd_cpu
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/retry_flat_same_dt.cpp -o build/retry_flat_same_dt
	build/retry_flat_same_dt
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/negative_pressure_policy.cpp -o build/negative_pressure_policy
	build/negative_pressure_policy
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/modern_dual_energy.cpp -o build/modern_dual_energy
	build/modern_dual_energy
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/parallel_shock.cpp -o build/parallel_shock
	build/parallel_shock
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/periodic_mean_field.cpp -o build/periodic_mean_field
	build/periodic_mean_field
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/brio_dual_energy.cpp -o build/brio_dual_energy
	build/brio_dual_energy
	python3 tests/check_run_sync.py
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/staggered_cfl.cpp -o build/staggered_cfl
	build/staggered_cfl
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/face_uct_1d.cpp -o build/face_uct_1d
	build/face_uct_1d
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/face_uct_2d.cpp -o build/face_uct_2d
	build/face_uct_2d
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/spherical_ballistic.cpp -o build/spherical_ballistic
	build/spherical_ballistic
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/boundary_contract.cpp -o build/boundary_contract
	build/boundary_contract
	$(CXX) $(CXXFLAGS) -I$(NATIVE)/include tests/transport_diagnostics.cpp -o build/transport_diagnostics
	build/transport_diagnostics
