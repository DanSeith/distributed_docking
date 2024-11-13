# Distributed Docking Pipeline

## Overview
This project aims to build a system that efficiently distributes docking jobs across workers. The input consists of a library of compounds (in SMILES/SDF format) and protein structures.

## Features
- **Queue Management**: Efficient handling of docking jobs.
- **Results Database**: Optimized storage and retrieval of results.
- **Status Tracking**: Monitoring job status and recovery from failures.
- **Smart Batching**: Grouping similar ligands for processing.
- **Caching**: Storing protein preparation steps to enhance performance.
- **API**: For job submission and monitoring.