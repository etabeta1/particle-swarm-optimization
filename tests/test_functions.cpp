#include "tests.hpp"
#include <point.hpp>
#include <function.hpp>

int main()
{
    TestManager m;

    m.addTest("SphereFunction at origin", []()
              {
        using namespace Swarm;
        std::vector<float> vals{0.0f, 0.0f};
        Point<float, 2> p(vals);
        Swarm::ObjectiveFunctions::SphereFunction<float, 2> f;
        return f.evaluate(p) == 0.0f ? PASS : ERR; });

    m.addTest("EllipsoidFunction at origin", []()
              {
        using namespace Swarm;
        std::vector<float> vals{0.0f, 0.0f};
        Point<float, 2> p(vals);
        ObjectiveFunctions::EllipsoidFunction<float, 2> f;
        return f.evaluate(p) == 0.0f ? PASS : ERR; });

    m.addTest("QuinticFunction at -1", []()
              {
        using namespace Swarm;
        std::vector<float> vals{-1.0f, -1.0f};
        Point<float, 2> p(vals);
        ObjectiveFunctions::QuinticFunction<float, 2> f;
        return f.evaluate(p) < 1e-5f ? PASS : ERR; });

    m.addTest("DropwaveFunction at origin", []()
              {
        using namespace Swarm;
        std::vector<float> vals{0.0f, 0.0f};
        Point<float, 2> p(vals);
        ObjectiveFunctions::DropwaveFunction<float, 2> f;
        return std::abs(f.evaluate(p)) < 1e-5f ? PASS : ERR; });

    m.addTest("Alpine1Function at origin", []()
              {
        using namespace Swarm;
        std::vector<float> vals{0.0f, 0.0f};
        Point<float, 2> p(vals);
        ObjectiveFunctions::Alpine1Function<float, 2> f;
        return std::abs(f.evaluate(p)) < 1e-5f ? PASS : ERR; });

    m.addTest("AckleyFunction at origin", []()
              {
        using namespace Swarm;
        std::vector<float> vals{0.0f, 0.0f};
        Point<float, 2> p(vals);
        ObjectiveFunctions::AckleyFunction<float, 2> f;
        return std::abs(f.evaluate(p)) < 1e-4f ? PASS : ERR; });

    m.addTest("SphereFunction at (1,2)", []()
              {
        using namespace Swarm;
        std::vector<float> vals{1.0f, 2.0f};
        Point<float, 2> p(vals);
        ObjectiveFunctions::SphereFunction<float, 2> f;
        return std::abs(f.evaluate(p) - 5.0f) < 1e-5f ? PASS : ERR; });

    m.addTest("EllipsoidFunction at (1,2)", []()
              {
        using namespace Swarm;
        std::vector<float> vals{1.0f, 2.0f};
        Point<float, 2> p(vals);
        ObjectiveFunctions::EllipsoidFunction<float, 2> f;
        return std::abs(f.evaluate(p) - 4.0f) < 1e-5f ? PASS : ERR; });

    if (m.runAll() == ERR)
        return 1;

    return 0;
}