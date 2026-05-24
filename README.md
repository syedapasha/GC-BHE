<div align="center">

# GC-BHE

### Uncover hidden causal links in event streams with Granger Causality testing for bivariate point processes.

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)](https://github.com/user/repo/actions)
[![License](https://img.shields.io/github/license/user/repo)](LICENSE)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](https://github.com/user/repo/pulls)
[![GitHub Stars](https://img.shields.io/github/stars/user/repo)](https://github.com/user/repo/stargazers)

</div>

---

## The Strategic "Why" (Overview)

> Traditional causality testing often falls short when dealing with dynamic, asynchronous event streams, leading to ambiguous insights and unreliable predictions in fields ranging from neuroscience to finance. Analyzing Granger Causality in bivariate point processes presents significant statistical and computational challenges, making it difficult for researchers and practitioners to accurately infer lead-lag relationships and predict future events based on past observations.

GC-BHE addresses these challenges by offering a specialized and robust R-based framework for rigorously testing Granger Causality within bivariate point processes. This project allows users to uncover causal relationships in complex event data, providing superior analytical depth and enabling more accurate predictive modeling and hypothesis testing than generic statistical methods.

---

## Key Features

✨ **Robust Granger Causality Inference**: Employ statistical methods tailored for point processes to provide high-confidence causal insights.
🔬 **Specialized Bivariate Point Process Analysis**: Designed from the ground up to handle the unique characteristics of event-based data, ensuring accurate and relevant results.
📊 **Integrated Data Generation & Simulation**: Includes utilities (`data_gen.r`) to create synthetic point process data, ideal for testing, validation, and understanding model behavior.
📈 **Expectation-Maximization (EM) Algorithm Implementation**: Utilizes EM algorithms (`em.r`) for robust parameter estimation in complex statistical models, enhancing model accuracy and convergence.
🎨 **Comprehensive Visualization Tools**: Generate insightful plots (`plots.r`, `gen_gc_plots.r`) to visually interpret Granger Causality test results and model diagnostics, making complex data understandable.
🔍 **Thorough Residual Analysis**: Provides tools (`residual.r`) for evaluating model fit and validating assumptions, ensuring the reliability of causal inferences.
🧩 **Modular R Script Design**: A clear and well-structured codebase (`main.r`, `test.r`) promotes ease of use, extensibility, and integration into existing R workflows.

---

## Technical Architecture

This project leverages the power of R for statistical computing, providing a robust environment for complex data analysis.


### Directory Structure

```
.
├── 📄 data_gen.r
├── 📄 em.r
├── 📄 gen_gc_plots.r
├── 📄 main.r
├── 📄 plots.r
├── 📄 residual.r
├── 📄 test.r
└── 📄 README.md
```

---

## Operational Setup

### Prerequisites

Ensure you have the R statistical programming language installed on your system.

*   **R**: [Download and Install R](https://cran.r-project.org/) (Version 4.0 or higher recommended)

### Installation

1.  **Clone the Repository**:
    ```bash
    git clone https://github.com/your-username/GC-BHE.git
    cd GC-BHE
    ```

2.  **Open R and Set Working Directory**:
    Launch RStudio or your preferred R environment, then set the working directory to the cloned repository:
    ```R
    # In R console
    setwd("/path/to/your/GC-BHE")
    ```

3.  **Install Required R Packages (if any)**:
    While specific packages are not listed in the manifest, R projects often depend on external libraries. If you encounter errors when running the scripts, you may need to install common data science packages. For instance:
    ```R
    # Example: Install commonly used packages (uncomment and run if needed)
    # install.packages(c("ggplot2", "dplyr", "tidyr", "lubridate"))
    ```

4.  **Run the Main Script or Tests**:
    Execute the primary script or the test suite to get started with the analysis:
    ```R
    # To run the main application logic or an example analysis
    source("main.r")

    # To run the test suite and validate functionality
    source("test.r")
    ```

---

