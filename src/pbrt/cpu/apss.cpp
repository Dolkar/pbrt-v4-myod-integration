#include <pbrt/cpu/apss.h>

#include <pbrt/shapes.h>
#include <pbrt/util/mesh.h>

namespace pbrt {

APSSPrimitive::APSSPrimitive(PointSet pointSet, Material material, int maxPrimsInNode,
                             BVHAggregate::SplitMethod splitMethod)
    : pointSet(pointSet), material(material), maxPrimsInNode(maxPrimsInNode) {
    primitiveMemory += sizeof(*this);
}

Bounds3f APSSPrimitive::Bounds() const {
    return Bounds3f();
}

pstd::optional<ShapeIntersection> APSSPrimitive::Intersect(const Ray& r,
                                                           Float tMax) const {
    return {};
}

bool APSSPrimitive::IntersectP(const Ray& r, Float tMax) const {
    return false;
}

}  // namespace pbrt