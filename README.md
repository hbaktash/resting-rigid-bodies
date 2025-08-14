# Putting Rigid Bodies to Rest

Analyzing and visualizing the resting behavior of rigid bodies under gravity with negligible momentum.

This repo focuses on:
- Computing the probability that an object rests on each stable face (via Morse-Smale complex of the gravitational potential function of a shape).
- Generating the quasi-static “drop” trajectory from a given initial orientation and exporting the rigid transforms.

Project page: https://hbaktash.github.io/projects/putting-rigid-bodies-to-rest/index.html

## Build

Requirements
- C++20 compiler (Clang 12+/GCC 10+)
- CMake ≥ 3.16
- Eigen3, Qhull
- (Optional) Polyscope

Install system packages
- macOS (Homebrew):
  - brew install cmake eigen qhull
- Ubuntu/Debian:
  - sudo apt-get update
  - sudo apt-get install -y build-essential cmake libeigen3-dev libqhull-dev

Clone
```
git clone https://github.com/hbaktash/Putting_Rigid_Bodies_to_Rest.git
cd Resting_Rigid_Bodies
```

Configure and build
- Headless (no visualization):
```
cmake -S . -B build -DRESTING_RIGID_BODIES_USE_POLYSCOPE=OFF
cmake --build build -j
```
- With visualization (Polyscope):
```
cmake -S . -B build -DRESTING_RIGID_BODIES_USE_POLYSCOPE=ON
cmake --build build -j
```

This produces the executable:
- build/drop_probs

CMake options
- RESTING_RIGID_BODIES_USE_POLYSCOPE=ON|OFF (default ON): enable/disable visualization features. When OFF, commands that request visualization are ignored with a warning.
- Inverse-design and other optional components are not required for the functionality documented here.

## Usage

The drop_probs tool supports:
- Probability mode: compute resting probability for each stable convex-hull face and write JSON.
- Drop mode: compute quasi-static transforms from an initial orientation and write JSON.
- Optional visualization: show the input orientation and the quasi-static trajectory on the Gauss map (requires a build with Polyscope).

Common flags
- --mesh <path>: input triangle mesh file (OBJ, etc.). The code orients faces, compresses the mesh, and uses the convex hull for analysis.
- --out <file>: output JSON file. Parent directories are created if needed.
- --viz: enable visualization for supported modes (requires building with Polyscope).
- Orientation (for drop mode):
  - --ox <float> --oy <float> --oz <float>: initial orientation direction (will be normalized). “Down” is (0, -1, 0).
- Mode selection:
  - --probs: probability mode
  - --drop: quasi-static drop mode

### Probability mode

Compute the probability mass associated with each stable convex-hull face. Probabilities are the face-region areas on the Gauss map divided by 4π.

Example
```
./build/drop_probs \
  --mesh data/dice.obj \
  --probs \
  --out results/dice_probabilities.json
```

Output (JSON)
```json
{
  "type": "stable_face_probabilities",
  "mesh_title": "dice",
  "stable_face_count": 6,
  "stable_faces": [
    { "face_index": 12, "normal": [0.0, -1.0, 0.0], "probability": 0.167 }
    // ...
  ]
}
```

### Drop mode

Generate the quasi-static “drop” sequence (rigid transforms) from an input orientation.

Example (headless)
```
./build/drop_probs \
  --mesh data/dice.obj \
  --ox 0 --oy -1 --oz 0 \
  --drop \
  --out results/dice_drop.json
```

Output (JSON)
```json
{
  "type": "quasi_static_drop",
  "mesh_title": "dice",
  "step_count": 42,
  "steps": [
    {
      "index": 0,
      "matrix": [
        [ 1.0, 0.0, 0.0, 0.0 ],
        [ 0.0, 1.0, 0.0, 0.0 ],
        [ 0.0, 0.0, 1.0, 0.0 ],
        [ 0.0, 0.0, 0.0, 1.0 ]
      ],
      "orientation": [0.0, -1.0, 0.0]
    }
    // ...
  ]
}
```

### Visualization (optional)

To visualize the input orientation and the quasi-static trajectory on the Gauss map:
1) Build with Polyscope:
```
cmake -S . -B build -DRESTING_RIGID_BODIES_USE_POLYSCOPE=ON
cmake --build build -j
```
2) Add --viz to your command:
```
./build/drop_probs \
  --mesh data/dice.obj \
  --ox 0 --oy -1 --oz 0 \
  --drop \
  --out results/dice_drop.json \
  --viz
```

Notes
- If built without Polyscope, --viz is ignored with a warning.
- Output files are created if they don’t exist; parent directories are also created.

## Data and formats

- Input: triangle mesh (e.g., OBJ). The pipeline orients faces, compresses, and uses its convex hull.
- Outputs:
  - Probability mode: JSON with per-face normals and probabilities.
  - Drop mode: JSON with ordered 4×4 transformation matrices and corresponding orientations.

## Citation

If you use this code in academic work, please cite the project as described on the project webpage:
- Project webpage: https://example.com/resting_rigid_bodies

## License

This code is intended for research and educational use. See