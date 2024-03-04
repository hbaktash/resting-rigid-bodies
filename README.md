# Rolling Dragons
Determining objects orientation when falling on the ground. Assuming no momentum.




### Get the code
Clone the project 
```
git clone --recursive git@git.corp.adobe.com:hbaktash/rolling-dragons-2D.git
```

#### Bullet library
For the bullet simulation executable I followed the instructions from https://github.com/erwincoumans/hello_bullet_cmake
Using a simpler build which does not handle MVS stuff.
Manually building Bullet is the way to go at the moment.

```
cd deps/bullet3
mkdir build_cmake
cd build_cmake
cmake -DCMAKE_INSTALL_PREFIX:PATH=local_install -DUSE_DOUBLE_PRECISION=ON -DCMAKE_DEBUG_POSTFIX="" -DINSTALL_LIBS=ON -DCMAKE_BUILD_TYPE=Release ..
make -j
make install
```

### Build the code

**Unix-like machines**: configure (with cmake) and compile
```
cd rolling-dragons-2D
mkdir build
cd build
cmake ..
make -j6
```

## Running the code

### polyhedra
This executable shows our implementation for the momentum-less case of dropping an object. All the visuals are included and can be worked with in Polyscope.

#### Note
* For building the raster image and building/splitin the Markov chain, center of mass needs to be inside the polyhedron; the check has been disabled for speed and can be enabled.
* Real-time generation for raster image can slow down the code depending on sample count. 
* Markov chain split can be done in real time without as much slow-down.

### bullet simulation
This is a simple simulation for dropping an objct on the ground. Momentum is not disabled yet and bounces, etc.. happen.
 

## More build notes

- manually pulled Polyscope for a recent update (initally detached at ..?)
- manually pulled imgui to be able to use https://github.com/epezent/implot/tree/master (initally detached at 512c54bb)