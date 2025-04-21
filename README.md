# P2Pool Sharechain Simulation Competency Test

## Overview
This repository contains the competency test implementation for the Summer of Bitcoin project proposal: "NS3 Simulation for P2Poolv2's Sharechain With Uncles." The simulation demonstrates a mesh/random P2P network with share distribution using NS-3.

## Implementation Details
The simulation creates a configurable P2P network with the following features:
- Configurable number of nodes (tested with 10+ nodes)
- Random peer connections with parameterized min/max peers per node
- "Share" message generation and propagation between peers
- Configurable share production intervals using normal distribution
- Network latency simulation with configurable parameters
- Comprehensive metrics collection on share propagation

## Video Demonstration
A video demonstration of the simulation using NS-3's NetAnim visualizer
[(demo)](https://github.com/user-attachments/assets/b9ad845a-199b-4a2b-8f32-22759f74dcb3)
The visualization shows:
- The network topology with connected peers
- Share message propagation across the network
- Node activity as shares are produced and relayed


## Running the Simulation
1. Ensure NS-3 is installed on your system
2. Clone this repository to your NS-3 workspace
3. Compile with `./waf --run "mesh-p2p-shares --nodes=12 --minPeers=2 --maxPeers=5"`
