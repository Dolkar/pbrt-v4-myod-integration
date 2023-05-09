#include <pbrt/cpu/apss.h>

#include <pbrt/shapes.h>
#include <pbrt/util/mesh.h>

namespace pbrt {


STAT_RATIO("Geometry/Points per set", nPts, nPtSets);

PointSet::PointSet(const Transform &renderFromObject, std::vector<Point3f> p,
                   std::vector<Normal3f> n, std::vector<Float> r, Allocator alloc)
    : nPoints(p.size()) {
    ++nPtSets;
    nPts += nPoints;

    // Transform points to rendering space and initialize point set _p_
    for (Point3f &pt : p)
        pt = renderFromObject(pt);
    this->p = point3BufferCache->LookupOrAdd(p, alloc);

    CHECK_EQ(nPoints, n.size());
    for (Normal3f &nn : n) {
        nn = renderFromObject(nn);
    }
    this->n = normal3BufferCache->LookupOrAdd(n, alloc);

    // Calculate uniform scale
    Float la2 = LengthSquared(renderFromObject(Vector3f(1, 0, 0)));
    Float lb2 = LengthSquared(renderFromObject(Vector3f(0, 1, 0)));
    Float lc2 = LengthSquared(renderFromObject(Vector3f(0, 0, 1)));
    Float scale = SafeSqrt((la2 + lb2 + lc2) / 3.0f);

    CHECK_EQ(nPoints, r.size());
    for (Float &rr : r)
        rr = scale * rr;
    this->r = floatBufferCache->LookupOrAdd(r, alloc);

    // Make sure that we don't have too much stuff to be using integers to
    // index into things.
    CHECK_LE(p.size(), std::numeric_limits<int>::max());
}

std::string PointSet::ToString() const {
    std::string np = "(nullptr)";
    return StringPrintf("[ PointSet nPoints: %d p: %s n: %s r: %s]", nPoints,
                        p ? StringPrintf("%s", pstd::MakeSpan(p, nPoints)) : np,
                        n ? StringPrintf("%s", pstd::MakeSpan(n, nPoints)) : np,
                        r ? StringPrintf("%s", pstd::MakeSpan(r, nPoints)) : np);
}

PointSet* PointSet::FromTriangleMesh(const TriangleMesh *mesh, float rMult,
                                    Allocator alloc) {
    // Assign a radius to each point as a function of the average distance to nearby points
    std::vector<Float> radii;
    radii.resize(mesh->nVertices);
    std::vector<int> degree;
    degree.resize(mesh->nVertices);

    for (int i = 0; i < mesh->nTriangles; i++) {
        // Update all three points
        int id0 = mesh->vertexIndices[i * 3 + 0];
        int id1 = mesh->vertexIndices[i * 3 + 1];
        int id2 = mesh->vertexIndices[i * 3 + 2];

        Point3f v0 = mesh->p[id0];
        Point3f v1 = mesh->p[id1];
        Point3f v2 = mesh->p[id2];

        Float d01 = Distance(v0, v1);
        Float d02 = Distance(v0, v2);
        Float d12 = Distance(v1, v2);

        radii[id0] += d01 + d02;
        radii[id1] += d01 + d12;
        radii[id2] += d02 + d12;

        degree[id0] += 2;
        degree[id1] += 2;
        degree[id2] += 2;
    }

    for (size_t i = 0; i < radii.size(); i++) {
        radii[i] = (radii[i] / degree[i]) * rMult;
    }

    std::vector<Point3f> p = std::vector<Point3f>(mesh->p, mesh->p + mesh->nVertices);
    std::vector<Normal3f> n = std::vector<Normal3f>(mesh->n, mesh->n + mesh->nVertices);

    return new PointSet({}, p, n, radii, alloc);
}

APSSPrimitive::APSSPrimitive(PointSet pointSet, Material material, int maxPrimsInNode,
                             BVHAggregate::SplitMethod splitMethod)
    : pointSet(pointSet), material(material) {
    primitiveMemory += sizeof(*this);

    // Build BVH from points
    std::vector<Primitive> inputPrimitives;
    for (int i = 0; i < pointSet.nPoints; i++) {
        inputPrimitives.push_back(new PointPrimitive(pointSet.p[i], pointSet.r[i], i));
    }

    BVHAggregate bvh = BVHAggregate(inputPrimitives, maxPrimsInNode, splitMethod);

    for (const Primitive &prim : bvh.GetPrimitives()) {
        primitives.push_back(*prim.Cast<PointPrimitive>());
    }
    nodes = bvh.GetFlattenedNodes();
}

Bounds3f APSSPrimitive::Bounds() const {
    CHECK(!nodes.empty());
    return nodes[0].bounds;
}

pstd::optional<ShapeIntersection> APSSPrimitive::Intersect(const Ray& r,
                                                           Float tMax) const {
    return {};
}

bool APSSPrimitive::IntersectP(const Ray& r, Float tMax) const {
    return false;
}

}  // namespace pbrt