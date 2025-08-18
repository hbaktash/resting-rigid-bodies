# Putting Rigid Bodies to Rest

Project page: https://hbaktash.github.io/projects/putting-rigid-bodies-to-rest/index.html

Analyzing and visualizing the resting behavior of rigid bodies under gravity with negligible momentum.

This repo focuses on:
- Computing the probability that an object rests on each stable orientation.
- Generating the quasi-static “drop” orientation trajectory from a given initial orientation.

For inverse design functionalities have a look at the [inverse design README](https://github.com/hbaktash/resting-rigid-bodies/blob/main/README_inverse_design).
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
git clone https://github.com/hbaktash/resting-rigid-bodies.git
cd resting-rigid-bodies
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
  "stable_face_count": 12,
  "stable_faces": [
    {
      "face_normal": [0.0,0.0,-1.0],
      "hull_face_index": 0,
      "probability": 0.0833333,
      "transformation_matrix": [
        [ 1.0, 0.0, 0.0, -0.0 ],
        [ 0.0, 0.0, 1.0, 0.5773],
        [ 0.0, -1.0, 0.0, -0.0],
        [ 0.0, 0.0, 0.0, 1.0]
        ]
    },
    //..
  ]
}
```

### Drop mode

Generate the quasi-static “drop” sequence (rigid transforms) from an input orientation. The ground is assumed to be at height 0.

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
  "final_orientation": [0.0, -1.0, 0.0],
  "initial_orientation": [0.0, -1.0, 0.0],
  "step_count": 1,
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
### Refining the Drop Sequence
By default the drop sequnce will include orientations where a 
new contact is made with the ground (e.g. vertex contact changes to edge contact). You can change this using the ```--min_QS_angle`` argument, 
which determines the maximum change in rotation angle at every step. 
By default this is set to 1 (radians), which wouldn't refine the rotation sequence.
Setting it to a lower value (0.001) will make the drop sequence look more smooth. 

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

## Citation

If this code contributes to academic work, please cite as:
```
@article{Baktash:2025:PRB,
  author      = {Baktash, Hossein and Sharp, Nicholas and
                 Zhou, Qingnan and Crane, Keenan and Jacobson, Alec},
  title       = {Putting Rigid Bodies to Rest},
  journal     = {ACM Trans. Graph.},
  issue_date  = {August&nbsp;2025},
  volume      = {44},
  number      = {4},
  articleno   = {155}, 
  numpages    = {16},
  year        = {2025},
  publisher   = {Association for Computing Machinery},
  address     = {New&nbsp;York, NY, USA},
  issn        = {0730-0301},
  url         = {https://doi.org/10.1145/3731203},
  doi         = {10.1145/3731203}
}
```
## License

This code is intended for research and educational use. See LICENSE for terms.