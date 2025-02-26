# ECG2O
**Equality-Constrained g2o – Extending g2o with Equality Constraint Support**

ECG2O extends the open-source **g2o (General Graph Optimization)** framework by introducing **equality constraints**. 
This allows for the addition of **equality edge factors**, enabling constrained optimization within factor graphs. 
The library also includes an **Augmented Lagrangian implementation** as a benchmark for constrained optimization techniques.

---

## Prerequisites

Before building and using ECG2O, ensure you have the following dependencies installed:

- **g2o:** Follow the official [g2o installation guide](https://github.com/RainerKuemmerle/g2o) to install the base library.
- **CMake:** Required for compiling the library. Install it via:
  ```bash
  sudo apt-get install cmake  # Ubuntu/Debian
  brew install cmake  # macOS
  ```

---

## Building ECG2O

1. **Navigate to the library folder:**
   ```bash
   cd path/to/ecg2o
   ```

2. **Create and enter a build directory:**
   ```bash
   mkdir build  
   cd build  
   ```

3. **Compile the library using CMake:**
   ```bash
   cmake ..  
   make  
   ```

---

## Usage

Once the library is successfully built, you can execute the provided functions for equality-constrained optimization.

```bash
./example_gn
```
or
```bash
./example_al
```
or
```bash
./example_soft
```
Ensure that all necessary input files or configurations are correctly set up before execution.

---

## License

This project is open-source and distributed under GPL-3.0
.
