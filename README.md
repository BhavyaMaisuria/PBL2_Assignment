# GenoMatch Pro 🧬

GenoMatch Pro is a high-performance C++ tool designed for genomic data analysis, focusing on DNA sequence alignment and fragment assembly. It features a modern, multi-threaded CLI with customizable UI themes.

## 🚀 Features

- **Global Alignment (Evolution Tracker):** Implements the Needleman-Wunsch algorithm to compare entire DNA sequences, ideal for tracking evolutionary changes between species.
- **Local Alignment (Mutation Scanner):** Uses the Smith-Waterman algorithm to identify highly similar local regions within DNA strands, perfect for mutation scanning and viral signature detection.
- **Fragment Assembly (Kahn's DAG Algorithm):** Reconstructs original DNA sequences from short fragments using a directed acyclic graph (DAG) and Kahn's topological sorting algorithm.
- **Dynamic UI Themes:** Switch between "BIO-CYAN" and "LAVA-RED" themes with a real-time pulsing UI effect powered by background threads.
- **Automated Data Loading:** Seamlessly syncs with `dna_data.txt` to load species and fragment data on startup.

## 🛠️ Algorithms Used

1.  **Needleman-Wunsch:** For global sequence alignment.
2.  **Smith-Waterman:** For local sequence alignment.
3.  **Kahn's Algorithm:** For topological sorting in fragment assembly.
4.  **Multi-threading:** Utilizes `std::thread` and `std::atomic` for background UI animations.

## 📂 Project Structure

- `atgc.cpp`: The core source code containing the logic for alignment, assembly, and the UI engine.
- `dna_data.txt`: The database containing genomic sequences and fragments for analysis.
- `atgc.exe`: Pre-compiled executable for Windows systems.

## 🚦 Getting Started

### Prerequisites

- A C++ compiler supporting C++11 or higher (e.g., GCC, Clang, or MSVC).

### Compilation

To compile the source code, use the following command:

```bash
g++ -o genomatch atgc.cpp -pthread
```

### Running the Application

After compilation, run the executable:

```bash
./genomatch
```

Ensure `dna_data.txt` is in the same directory as the executable for proper data synchronization.

## 📊 Data Format (`dna_data.txt`)

The system reads data in a simplified format:
- Lines starting with `#` or `>` are treated as comments/headers and skipped.
- Long strings (>10 characters) are treated as full DNA strands for alignment.
- Short strings are treated as fragments for assembly.

## 🎨 UI Controls

- **1-3:** Select analysis algorithms.
- **4:** Toggle between BIO-CYAN and LAVA-RED themes.
- **0:** Safely shutdown the system.
