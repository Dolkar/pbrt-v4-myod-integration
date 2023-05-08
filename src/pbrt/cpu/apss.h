#ifndef PBRT_CPU_APSS_H
#define PBRT_CPU_APSS_H

#include <pbrt/pbrt.h>

#include <pbrt/cpu/aggregates.h>
#include <pbrt/util/mesh.h>

namespace pbrt {

// #MYOD
// APSSPrimitive Definition
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
    int maxPrimsInNode;
};

}  // namespace pbrt

#endif  // PBRT_CPU_APSS_H