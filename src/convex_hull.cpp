/************************************************************************
*
* ADOBE CONFIDENTIAL
* ___________________
*
* Copyright [first year code created] Adobe
* All Rights Reserved.
*
* NOTICE: All information contained herein is, and remains
* the property of Adobe and its suppliers, if any. The intellectual
* and technical concepts contained herein are proprietary to Adobe
* and its suppliers and are protected by all applicable intellectual
* property laws, including trade secret and copyright laws.
* Dissemination of this information or reproduction of this material
* is strictly forbidden unless prior written permission is obtained
* from Adobe.
*************************************************************************
*/
#include "convex_hull.h"

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> 
convex_hull(std::vector<Vector3> point_set){
    // for (Vector3 p: point_set){
    //     std::vector<coordT> p_vec({p.x, p.y, p.z});
    //     rbox.append(orgQhull::PointCoordinates(p_vec, 3));
    // }
    const size_t num_points = point_set.size();
    const size_t dim = 3;
    char flags[64];
    sprintf(flags, "qhull Qt");
    coordT* data = new coordT[num_points * dim];
    std::copy(point_set.data(), point_set.data() + num_points * dim, data);
    RboxPoints rbox;
    rbox.append(orgQhull::PointCoordinates()
    Qhull qhull;
    
    
    // int err = qh_new_qhull(dim, num_points, data, false,
    //         flags, NULL, stderr);

    // if (!err) {
    //     extract_hull(points);
    // } else {
    //     std::stringstream err_msg;
    //     err_msg << "Qhull error: " << err;
    //     throw RuntimeError(err_msg.str());
    // }

    // qh_freeqhull(!qh_ALL);
    // m_lock.unlock();
    // delete [] data;
    // reorient_faces();
}

void QhullEngine::extract_hull(const MatrixFr& points) {
    const size_t dim = points.cols();
    const size_t num_input_points = points.rows();
    const size_t num_faces = qh num_facets;
    const size_t num_vertices = qh num_vertices;

    size_t index = 0;
    m_vertices.resize(num_vertices, dim);
    m_index_map.resize(num_vertices);
    VectorF inverse_map = VectorF::Ones(num_input_points) * -1;
    vertexT* vertex, **vertexp;
    FORALLvertices {
        size_t i = qh_pointid(vertex->point);
        m_vertices.row(index) = Eigen::Map<VectorF>(vertex->point, dim);
        m_index_map[index] = i;
        inverse_map[i] = index;
        index++;
    }

    index = 0;
    m_faces.resize(num_faces, dim);
    facetT *facet;
    FORALLfacets {
        size_t i=0;
        FOREACHvertex_( facet->vertices ) {
            assert(inverse_map(qh_pointid(vertex->point)) != -1);
            m_faces(index, i) = inverse_map(qh_pointid(vertex->point));
            i++;
        }
        index++;
    }
}