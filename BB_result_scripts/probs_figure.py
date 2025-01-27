#!/usr/bin/env python

import hakowan as hkw
import lagrange
import numpy as np

title = "piggy"
mesh = lagrange.io.load_mesh(f"{title}.obj")
hull = lagrange.io.load_mesh(f"hull_{title}.obj")
lagrange.normalize_mesh(mesh)
lagrange.normalize_mesh(hull)
nid = lagrange.compute_facet_normal(hull)
hull_normals = hull.attribute(nid).data
#center_of_mass = np.array([-0.0137228, -0.01872048, -0.03269057])
center_of_mass = np.array([0, 0, 0])


def rotation(from_vector, to_vector):
    axis = np.cross(from_vector, to_vector)
    sin_a = np.linalg.norm(axis)
    cos_a = np.dot(from_vector, to_vector)
    if sin_a < 1e-9:
        return np.eye(4)
    else:
        v = np.array(axis / sin_a, dtype=np.float64)
        I = np.eye(3)
        H = np.outer(v, v)
        S = np.cross(I, v)
        M = I * cos_a + S * sin_a + H * (1 - cos_a)

        A = np.eye(4, dtype=np.float64)
        A[:3, :3] = M
        return A


def get_transformation(fid):
    n = hull_normals[fid]
    minus_y = np.array([0, -1, 0])

    # Move center of mass to origin
    t1 = np.eye(4)
    t1[:3, 3] = -center_of_mass

    # Rotate such that facet indexed by fid is facing the -y direction
    r = rotation(n, minus_y)

    # Move the shape up so that it rests on the floor
    d = np.dot(hull.vertices[hull.facets[fid][0]] - center_of_mass, n)
    t2 = np.eye(4)
    t2[:3, 3] = -d * minus_y

    return t2 @ r @ t1


def generate_floor(w, h, z_up):
    floor = lagrange.SurfaceMesh()
    z = 0
    vertices = np.array(
        [
            [-w / 2, -h / 2, z],
            [w / 2, -h / 2, z],
            [w / 2, h / 2, z],
            [-w / 2, h / 2, z],
        ]
    )
    if not z_up:
        vertices = vertices[:, [0, 2, 1]]
        vertices[:, 2] *= -1

    floor.add_vertices(vertices)
    floor.add_quad(0, 1, 2, 3)
    uvs = np.array(
        [
            [0, 0],
            [1, 0],
            [1, 1],
            [0, 1],
        ],
        dtype=np.float32,
    )
    floor.create_attribute(
        "uv",
        element=lagrange.AttributeElement.Vertex,
        usage=lagrange.AttributeUsage.UV,
        initial_values=uvs,
    )
    return floor


M_330 = get_transformation(330)   # prob 0.2078
M_130 = get_transformation(130)   # 0.2697
M_653 = get_transformation(653)   # 0.10195
M_31 = get_transformation(31)     # 0.31884
M_43 = get_transformation(43)     # 0.05448
M_1590 = get_transformation(1590) # 0.0381
M_3303 = get_transformation(3303) # 0.0.0073

object_layer = hkw.layer().data(mesh).material("RoughPlastic", "#3B95B3")

config_330 = object_layer.transform(hkw.transform.Affine(M_330))
config_130 = object_layer.transform(hkw.transform.Affine(M_130))
config_653 = object_layer.transform(hkw.transform.Affine(M_653))
config_31 = object_layer.transform(hkw.transform.Affine(M_31))
config_43 = object_layer.transform(hkw.transform.Affine(M_43))
config_1590 = object_layer.transform(hkw.transform.Affine(M_1590))
config_3303 = object_layer.transform(hkw.transform.Affine(M_3303))

frames = (
    config_330.translate([-3, 0, 0])
    + config_130.translate([-1, 0, 0])
    + config_653.translate([1, 0, 0])
    + config_31.translate([3, 0, 0])
    + config_43.translate([5, 0, 0])
    + config_1590.translate([7, 0, 0])
    + config_3303.translate([9, 0, 0])
)

roi = np.array(
    [
        [-3, -0.5, -1],
        [3, 1.5, 1],
    ],
    dtype=np.float32,
)

floor = generate_floor(20, 20, False)
floor_layer = (
    hkw.layer()
    .data(floor, roi_box=roi)
    .material(
        "Principled",
        # "lightgray",
        hkw.texture.Checkerboard(
            texture1=0.5,
            texture2=0.8,
            size=8,
        ),
        two_sided=True,
    )
)

config = hkw.config()
#config.z_up()
config.sensor.location = [0, 0, 2]
config.film.width = 1920
config.film.height = 800

hkw.render(
   frames + floor_layer,
   config,
   filename=f"{title}_probs_figure.png",
)