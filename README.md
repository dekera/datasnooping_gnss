# Geodetic Adjustment with Outlier Detection (IDS)

This repository contains a Python implementation for performing geodetic network adjustment using the Least Squares method, with an added Iterative Data Snooping (IDS) algorithm for outlier detection and removal. This method is designed to identify and remove anomalous observations in geodetic networks, improving the accuracy and reliability of the adjusted coordinates.

## Overview

This project focuses on adjusting a geodetic network composed of reference points and baselines between them. The iterative process detects outliers based on residuals, and once identified, the outliers are removed to enhance the precision of the adjustment. The outlier detection is performed using the **Iterative Data Snooping (IDS)** algorithm, which iteratively removes the observation with the largest residual until no more outliers are detected or the maximum allowable residual is below a defined threshold.

The adjustment method is based on **Least Squares (LS)**, a statistical technique used to minimize the sum of squared residuals between the observed and calculated measurements. The IDS method identifies outliers by calculating the normalized residuals and removing data points that significantly deviate from the expected values.

### Key Features:

* Geodetic adjustment using the **Least Squares** method.
* Iterative outlier detection using the **IDS** algorithm.
* Generation of a comprehensive report detailing the removed outliers and adjusted coordinates.
* Capability to work with 3D geodetic coordinates and multiple baselines between points.

## Prerequisites

This project uses **Pixi** to manage dependencies. Ensure that **Pixi** is installed to handle the project’s requirements.

To install **Pixi**, follow the instructions on the [Pixi documentation](https://pixi.com).

### Installation and Setup

1. **Clone the repository**:

```bash
git clone https://github.com/YOUR_USERNAME/YOUR_REPOSITORY.git
cd YOUR_REPOSITORY
```

2. **Install dependencies** using Pixi:

```bash
pixi install
```

3. **Run the script** to perform geodetic adjustment and outlier detection:

```bash
python datasnooping.py
```

The script will execute the geodetic adjustment, apply IDS to detect and remove outliers, and generate a report detailing the results.

## Usage

### Code Walkthrough

1. **Fixed Coordinates**: The coordinates of the reference point, BEPA, are predefined as fixed values. These values are used in the adjustment process as the baseline for calculating relative positions of other points in the network.

2. **Baseline Data**: The baseline data defines the distances between pairs of points in the geodetic network, along with their measurement uncertainties (sigma values). The baseline data includes both horizontal (X, Y) and vertical (Z) displacements.

3. **Matrix Construction**:

   * The matrix **A** (design matrix) is constructed based on the relationships between the points and the observations.
   * **L** (observation vector) contains the actual measurements, i.e., the differences in coordinates between the points.
   * **Sigma** represents the uncertainties associated with the measurements.

4. **Adjustment and IDS**:

   * The geodetic adjustment is performed using the Least Squares method, where the solution is obtained by minimizing the sum of squared residuals.
   * The IDS algorithm iteratively removes outliers based on the largest residuals. If the absolute residual exceeds a predefined critical value (3.29), the corresponding observation is flagged as an outlier and removed.

5. **Results**: After all iterations, the script outputs a final adjusted coordinate solution and a summary of the outliers removed during the process.

### Functions

* **coord_bepa(comp)**: Returns the fixed coordinates of the BEPA point (X, Y, Z).
* **add_eq(i, j, d, sigma_m, comp)**: Adds the equations corresponding to the baselines between two points to the design matrix **A**.
* **ajusta(A_sub, L_sub, sigma_sub)**: Performs the geodetic adjustment using the Least Squares method, computes residuals, and applies the IDS algorithm to detect outliers.

### Error Handling

The algorithm will print the **maximum residual** for each iteration and flag the corresponding outliers. The script will continue iterating until no further outliers are detected or the maximum residual is below the critical threshold.

## Dependencies

The following dependencies are required for this project:

* **Pixi**: Used for dependency management.
* **NumPy**: Used for matrix and array manipulations.

To install the dependencies with Pixi:

```bash
pixi install
```

## Conclusion

This repository provides a Python implementation for geodetic network adjustment with outlier detection. It is useful for geospatial applications that require precise network adjustments and the identification of anomalous measurements. The IDS algorithm ensures that the final adjusted coordinates are as accurate as possible by removing outliers that may distort the results.
