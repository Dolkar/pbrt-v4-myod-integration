#ifndef PBRT_CPU_APSS_H
#define PBRT_CPU_APSS_H

#include <pbrt/pbrt.h>

#include <pbrt/cpu/aggregates.h>
#include <pbrt/util/mesh.h>

namespace pbrt {

// #MYOD
struct PointSet {
  public:
    // PointSet Public Methods
    PointSet(const Transform &renderFromObject, std::vector<Point3f> p,
             std::vector<Normal3f> n, std::vector<Float> r, Allocator alloc);

    std::string ToString() const;

    static PointSet* FromTriangleMesh(const TriangleMesh* mesh, float rMult, Allocator alloc);

    // PointSet Public Members
    int nPoints;
    const Point3f *p = nullptr;
    const Normal3f *n = nullptr;
    const Float *r = nullptr;
};

// #MYOD
class PointPrimitive {
  public:
    // PointPrimitive Public Methods
    PointPrimitive(Point3f p, Float r, int id)
        : p(p), r(r), id(id) {}

    Bounds3f Bounds() const {
        return Bounds3f(p - Vector3f(r, r, r), p + Vector3f(r, r, r));
    }
    pstd::optional<ShapeIntersection> Intersect(const Ray &r, Float tMax) const { return {}; }
    bool IntersectP(const Ray &r, Float tMax) const { return false; }

  private:
    // PointPrimitive Private Members
    Point3f p;
    Float r;
    int id;
};

// #MYOD
class APSSPrimitive {
  public:
    // APSSPrimitive Public Methods
    APSSPrimitive(PointSet pointSet, Material material, int maxPrimsInNode = 1,
                  BVHAggregate::SplitMethod splitMethod = BVHAggregate::SplitMethod::SAH);
    Bounds3f Bounds() const;
    pstd::optional<ShapeIntersection> Intersect(const Ray &r, Float tMax) const;
    bool IntersectP(const Ray &r, Float tMax) const;

  private:
    // APSSPrimitive Private Members
    PointSet pointSet;
    Material material;
    std::vector<PointPrimitive> primitives;
    std::vector<LinearBVHNode> nodes;
};

}  // namespace pbrt

#endif  // PBRT_CPU_APSS_H