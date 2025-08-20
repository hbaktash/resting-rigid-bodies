# Putting Rigid Bodies to Rest

Project page: https://hbaktash.github.io/projects/putting-rigid-bodies-to-rest/index.html

Analyzing and visualizing the resting behavior of rigid bodies under gravity with negligible momentum.

This repo focuses on:
- Computing the probability that an object rests on each stable orientation.
- Generating the quasi-static “drop” orientation trajectory from a given initial orientation.

For inverse design functionalities have a look at the [inverse design README](https://github.com/hbaktash/resting-rigid-bodies/blob/main/README_inverse_design).

### Quick Examples

GUI example:
```bash
./build/drop_probs \
  --mesh data/tet.obj \
  --viz \
  --out ./tet_results.json
```

This will load the Polyscope GUI and provide the interface to initialize with an orientation and drop the mesh on the ground. The drop sequence can be visualized using the provided slider.

<p align="center">
  <img src="./data/drop_prob_UI.gif" width="500" alt="convex dice demo" />
</p>

The probabilities of each stable face, that correspond to red points on the Gauss Sphere (local minima of the potential function), can be seen by clicking on each read point.

You can skip the GUI and directly output the results to file using `--drop` and `--probs` flags. 

Saving probabilities to file:
```bash
./build/drop_probs \
  --mesh data/tet.obj \
  --probs \
  --out tet_results.json
```

Saving a drop sequence to file (need to provide initial orientation in command line):
```bash
./build/drop_probs \
  --mesh data/tet.obj \
  --drop \
  --ox 0 \
  --oy -1 \
  --oz 0 \
  --out tet_results.json
```

Saving probabilities of each stable orientation to file:
```bash
./build/drop_probs \
  --mesh data/tet.obj \
  --probs \
  --out ./tet_result.json
```


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
- With visualization (Polyscope):
```
cmake -S . -B build -DRESTING_RIGID_BODIES_USE_POLYSCOPE=ON
cmake --build build -j
```
- No visualization:
```
cmake -S . -B build -DRESTING_RIGID_BODIES_USE_POLYSCOPE=OFF
cmake --build build -j
```

This produces the executable:
- build/drop_probs

Other CMake options:
- Inverse-design functionality is documented in [README_inverse_design](https://github.com/hbaktash/resting-rigid-bodies/blob/main/README_inverse_design.md) .

## Usage

The drop_probs tool supports:
- Probability mode: compute resting probability for each stable convex-hull face.
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
- --min_QS_angle <float>: minimum angle change (radians) for each quasi-static step (default 1.0). Lower values yield smoother trajectories. The default is 1.0 and it only updates the orientation when a new contact is made with the ground.

### Probability mode

Compute the probability mass associated with each stable convex-hull face. Probabilities are the Morse-Smale complex cell areas on the Gauss map divided by 4π.

Example
```
./build/drop_probs \
  --mesh data/dice.obj \
  --probs \
  --out ./dice_result.json
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
  ]
}
```

The output includes:
- stable_face_count: number of stable faces (local minima on Gauss map)
- stable_faces: array of stable face data
  - face_normal: normal vector of the stable face
  - hull_face_index: index of the face in the convex hull mesh
  - probability: resting probability on this face
  - transformation_matrix: 4x4 matrix to align the face normal with the downward direction.

### Drop mode

Generate the quasi-static “drop” sequence (rigid transforms) from an input orientation. The ground is assumed to be at height 0.

Example (headless)
```
./build/drop_probs \
  --mesh data/dice.obj \
  --ox 0 --oy -1 --oz 0 \
  --drop \
  --out ./dice_result.json
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
    },
  ]
}
```
The output includes:
- final_orientation: final resting orientation (down vector)
- initial_orientation: input orientation (down vector)
- step_count: number of steps in the drop sequence
- steps: array of step data
  - index: step index
  - matrix: 4x4 transformation matrix for this step
  - orientation: orientation at this step (with respect to the initial reference pose)

### Refining the Drop Sequence
By default the drop sequence will include orientations where a
new contact is made with the ground (e.g. vertex contact changes to edge contact). You can change this using the ```--min_QS_angle``` argument,
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
  --out ./dice_result.json \
  --viz
```

An example usage of the GUI is shown above (quick example).

Notes
- If built without Polyscope, --viz is ignored with a warning.
- Output parent directories are created if they don’t exist.


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
This code is licensed under the MIT License (see LICENSE for details).