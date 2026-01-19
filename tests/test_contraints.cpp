#include "tests.hpp"
#include <constraint.hpp>

int main()
{
    TestManager m;

    m.addTest("Test box constraint inside", []()
              {
        using namespace Swarm;

        Point<float, 2> a(0.0f);
        Point<float, 2> b(10.0f);
        auto box_constraint = makeBoxConstraint(a, b);

        Point<float, 2> inside_point(5.0f);
        if (!box_constraint(inside_point))
            return ERR;
        return PASS; });

    m.addTest("Test box constraint outside", []()
              {
        using namespace Swarm;

        Point<float, 2> a(0.0f);
        Point<float, 2> b(10.0f);
        auto box_constraint = makeBoxConstraint(a, b);

        Point<float, 2> outside_point(15.0f);
        if (box_constraint(outside_point))
            return ERR;
        return PASS; });

    m.addTest("Test negate constraint", []()
              {
        using namespace Swarm;
        Point<float, 2> a(0.0f);
        Point<float, 2> b(10.0f);
        auto box_constraint = makeBoxConstraint(a, b);
        auto negated_constraint = negateConstraint(box_constraint);
        Point<float, 2> outside_point(15.0f);
        if (!negated_constraint(outside_point))
            return ERR;
        Point<float, 2> inside_point(5.0f);
        if (negated_constraint(inside_point))
            return ERR;
        return PASS; });

    m.addTest("Test max distance constraint", []()
              {
        using namespace Swarm;
        Point<float, 2> center(0.0f);
        float max_distance = 5.0f;
        auto max_distance_constraint = makeMaxDistanceConstraint(center, max_distance);

        std::vector<float> coords_inside = {3.0f, 4.0f};
        Point<float, 2> inside_point(coords_inside); // Distance is 5.0
        if (!max_distance_constraint(inside_point))
            return ERR;

        std::vector<float> coords_outside = {6.0f, 8.0f};
        Point<float, 2> outside_point(coords_outside); // Distance is 10.0
        if (max_distance_constraint(outside_point))
            return ERR;

        return PASS; });

    m.addTest("Test elliptic constraint", []()
              {
        using namespace Swarm;
        std::vector<float> f1_coords = {0.0f, 0.0f};
        std::vector<float> f2_coords = {4.0f, 0.0f};
        Point<float, 2> f1(f1_coords);
        Point<float, 2> f2(f2_coords);
        float tol = 6.0f;
        auto elliptic_constraint = makeEllipticConstraint(f1, f2, tol);

        std::vector<float> coords_inside = {2.0f, 1.0f};
        Point<float, 2> inside_point(coords_inside); // Sum of distances is approx 6.47
        if (!elliptic_constraint(inside_point))
            return ERR;

        std::vector<float> coords_outside = {5.0f, 5.0f};
        Point<float, 2> outside_point(coords_outside); // Sum of distances is approx 8.94
        if (elliptic_constraint(outside_point))
            return ERR;

        return PASS; });

    if (m.runAll() == ERR)
        return 1;

    return 0;
}