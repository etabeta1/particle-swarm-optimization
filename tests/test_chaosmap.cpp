#include "tests.hpp"
#include <chaosmap.hpp>

int main()
{
    TestManager m;

    m.addTest("Test custom ChaosMap", []()
              {
        using namespace Swarm;

        Point<float, 2> m_a(0.f);
        Point<float, 2> m_b(1.f);

        ChaosMap<float, 2> custom_map(
            [](const Point<float, 2> &p, int)
            {
                std::vector<float> coords = {0.1f, 0.2f};
                return p + Point<float, 2>(coords);
            }, m_a, m_b);

        Point<float, 2> p(7.5f);
        Point<float, 2> a(5.f);
        Point<float, 2> b(10.f);
        Point<float, 2> p_mapped = custom_map.getPoint(p, a, b, 1);

        if (std::abs(p_mapped[0] - 8.f) > 1e-4f || std::abs(p_mapped[1] - 8.5f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test ChaosFactory::Chebyshev", []()
              {
        using namespace Swarm;

        Point<float, 1> m_a(-1.0);
        Point<float, 1> m_b(1.0);

        ChaosFactory::Chebyshev<float, 1> chebyshev_map;

        Point<float, 1> p(1.0);
        Point<float, 1> a(-2.0);
        Point<float, 1> b(2.0);
        Point<float, 1> p_mapped = chebyshev_map.getPoint(p, a, b, 2);

        if (std::abs(p_mapped[0] + 1.f) > 1e-4)
            return ERR;
        return PASS; });

    m.addTest("Test ChaosFactory::Singer", []()
              {
        using namespace Swarm;

        Point<float, 1> m_a(0.0);
        Point<float, 1> m_b(1.0);

        ChaosFactory::Singer<float, 1> singer_map(0.9);

        Point<float, 1> p(0.5);
        Point<float, 1> a(0.0);
        Point<float, 1> b(2.0);
        Point<float, 1> p_mapped = singer_map.getPoint(p, a, b, 0);

        if (std::abs(p_mapped[0] - 1.62968) > 1e-4)
            return ERR;
        return PASS; });

    m.addTest("Test ChaosFactory::Sine", []()
              {
        using namespace Swarm;

        Point<float, 1> m_a(0.0);
        Point<float, 1> m_b(1.0);
        ChaosFactory::Sine<float, 1> sine_map(0.8);
        
        Point<float, 1> p(0.5);
        Point<float, 1> a(0.0);
        Point<float, 1> b(1.0);
        Point<float, 1> p_mapped = sine_map.getPoint(p, a, b, 0);
        const float pi = std::numbers::pi;

        float expected = 0.8 * std::sin(pi * 0.5);

        if (std::abs(p_mapped[0] - expected) > 1e-4)
            return ERR;
        return PASS; });

    m.addTest("Test ChaosFactory::Sinusoidal", []()
              {
        using namespace Swarm;

        Point<float, 1> m_a(0.0);
        Point<float, 1> m_b(1.0);
        ChaosFactory::Sinusoidal<float, 1> sinusoidal_map(0.9);
        
        Point<float, 1> p(0.5);
        Point<float, 1> a(0.0);
        Point<float, 1> b(1.0);
        Point<float, 1> p_mapped = sinusoidal_map.getPoint(p, a, b, 0);
        const float pi = std::numbers::pi;

        float expected = 0.9 * 0.5 * 0.5 * std::sin(pi * 0.5);

        if (std::abs(p_mapped[0] - expected) > 1e-4)
            return ERR;
        return PASS; });

    m.addTest("Test ChaosFactory::LogisticMap", []()
              {
        using namespace Swarm;
        Point<float, 1> m_a(0.0);
        Point<float, 1> m_b(1.0);
        ChaosFactory::LogisticMap<float, 1> logistic_map(3.5);

        Point<float, 1> p(0.5);
        Point<float, 1> a(0.0);
        Point<float, 1> b(1.0);
        Point<float, 1> p_mapped = logistic_map.getPoint(p, a, b, 0);

        float expected = 3.5 * 0.5 * (1 - 0.5);
        if (std::abs(p_mapped[0] - expected) > 1e-4)
            return ERR;
        return PASS; });

    m.addTest("Test ChaosFactory::Iterative", []()
              {
        using namespace Swarm;
        Point<float, 1> m_a(-1.0);
        Point<float, 1> m_b(1.0);
        ChaosFactory::Iterative<float, 1> iterative_map(2.0);

        Point<float, 1> p(0.5);
        Point<float, 1> a(-3.0);
        Point<float, 1> b(1.0);
        Point<float, 1> p_mapped = iterative_map.getPoint(p, a, b, 3);

        if (std::abs(p_mapped[0] - 0.73205) > 1e-4)
            return ERR;
        return PASS; });

    if (m.runAll() == ERR)
        return 1;

    return 0;
}