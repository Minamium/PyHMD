# Development contracts

- Keep the supported mathematical scope explicit: 2D Cartesian / spherical axisymmetric resistive MHD, FP64, existing settings only until a task requests an extension.
- Python owns configuration, launch and results; native C++/CUDA owns MPI, buffers, time integration and accepted-step accounting.
- Register executable capabilities only. Reject unsupported options instead of silently selecting a different backend or model.
- Separate numerical changes from structural changes. Preserve the existing CPU regression suite and add boundary/API tests when changing those contracts.
- Run `make test` and `make test-python` for relevant native/API changes. Verify packaging if bundled sources or build code change. Do not claim CUDA/MPI validation from a CPU run.
- Keep outputs and metadata honest: completed integration, numerical integrity, and scientific reproduction are distinct claims.
- Keep generated runs, binaries, credentials and personal cluster connection details out of commits.
- Preserve old result readability and strict checkpoint signature checks; do not silently bypass a generation mismatch.
