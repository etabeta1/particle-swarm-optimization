#include "tests.hpp"
#include <point.hpp>

int main()
{
    TestManager m;

    m.addTest("Test point vector constructor", []()
              {
        using namespace Swarm;

        std::vector<float> vec = {1.0f, 2.0f, 3.0f};
        Point<float, 3> p1(vec);
        if (p1[0] != 1.0f || p1[1] != 2.0f || p1[2] != 3.0f)
            return ERR;
        return PASS; });

    m.addTest("Test point constructors with array and scalar", []()
              {
        using namespace Swarm;

        float arr[3] = {4.0f, 5.0f, 6.0f};
        Point<float, 3> p2(arr);
        if (p2[0] != 4.0f || p2[1] != 5.0f || p2[2] != 6.0f)
            return ERR;
        return PASS; });

    m.addTest("Test point constructors with scalar", []()
              {
        using namespace Swarm;

        Point<float, 3> p3(7.0f);
        if (p3[0] != 7.0f || p3[1] != 7.0f || p3[2] != 7.0f)
            return ERR;
        return PASS; });

    m.addTest("Test point copy constructor", []()
              {
        using namespace Swarm;

        std::vector<float> vec = {1.0f, 2.0f, 3.0f};
        Point<float, 3> p1(vec);
        Point<float, 3> p4(p1);
        if (p4[0] != 1.0f || p4[1] != 2.0f || p4[2] != 3.0f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::operator+", []()
              {
        using namespace Swarm;

        Point<float, 2> p1(1.0f);
        Point<float, 2> p2(2.0f);
        Point<float, 2> p3 = p1 + p2;
        if (std::abs(p3[0] - 3.0f) > 1e-4f || std::abs(p3[1] - 3.0f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::operator-", []()
              {
        using namespace Swarm;

        Point<float, 2> p1(1.0f);
        Point<float, 2> p2(2.0f);
        Point<float, 2> p4 = p2 - p1;
        if (std::abs(p4[0] - 1.0f) > 1e-4f || std::abs(p4[1] - 1.0f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::operator* (with scalar)", []()
              {
        using namespace Swarm;

        Point<float, 2> p1(1.0f);
        Point<float, 2> p5 = p1 * 3.0f;
        if (std::abs(p5[0] - 3.0f) > 1e-4f || std::abs(p5[1] - 3.0f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::operator* (with other Point)", []()
              {
        using namespace Swarm;

        Point<float, 2> p1([](size_t index)
                           { return static_cast<float>(index + 1); });
        Point<float, 2> p2([](size_t index)
                           { return static_cast<float>(index + 2); });
        Point<float, 2> p3 = p1 * p2;
        if (std::abs(p3[0] - 2.0f) > 1e-4f || std::abs(p3[1] - 6.0f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::operator/ (with scalar)", []()
              {
        using namespace Swarm;

        Point<float, 2> p1(4.0f);
        Point<float, 2> p2 = p1 / 2.0f;
        if (std::abs(p2[0] - 2.0f) > 1e-4f || std::abs(p2[1] - 2.0f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::operator/ (with other Point)", []()
              {
        using namespace Swarm;

        Point<float, 2> p1([](size_t index)
                           { return static_cast<float>(index + 2); });
        Point<float, 2> p2([](size_t index)
                           { return static_cast<float>(index + 4); });
        Point<float, 2> p3 = p2 / p1;
        if (std::abs(p3[0] - 2.0f) > 1e-4f || std::abs(p3[1] - (5.f/3.f)) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::clamp", []()
              {
        using namespace Swarm;

        Point<float, 2> p1(-1.0f);
        Point<float, 2> p2(+3.0f);
        Point<float, 2> p3([](size_t index)
                           { return static_cast<float>(index*10) - 5; });
        Point<float, 2> p4 = p3.clamp(p1, p2);
        if (std::abs(p4[0] + 1.0f) > 1e-4f || std::abs(p4[1] - 3.0f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::pow", []()
              {
        using namespace Swarm;
        
        Point<float, 2> p5(3.0f);
        Point<float, 2> p6 = p5.pow(2);
        if (std::abs(p6[0] - 9.0f) > 1e-4f || std::abs(p6[1] - 9.0f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::sin", []()
              {
        using namespace Swarm;

        Point<float, 2> p1([](size_t index)
                           { return static_cast<float>(index * 3.14159265 / 2); });
        Point<float, 2> p2 = p1.sin();
        if (std::abs(p2[0] - 0.0f) > 1e-4f || std::abs(p2[1] - 1.0f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::arcsin", []()
              {
        using namespace Swarm;
        Point<float, 2> p1([](size_t index)
                           { return static_cast<float>(index); });
        Point<float, 2> p2 = p1.arcsin();
        if (std::abs(p2[0] - 0.0f) > 1e-4f || std::abs(p2[1] - 1.57079633f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::cos", []()
              {
        using namespace Swarm;
        Point<float, 2> p1([](size_t index)
                           { return static_cast<float>(index * 3.14159265 / 2); });
        Point<float, 2> p2 = p1.cos();
        if (std::abs(p2[0] - 1.0f) > 1e-4f || std::abs(p2[1] - 0.0f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::arccos", []()
              {
        using namespace Swarm;
        Point<float, 2> p1([](size_t index)
                           { return static_cast<float>(index); });
        Point<float, 2> p2 = p1.arccos();
        if (std::abs(p2[0] - 1.57079633f) > 1e-4f || std::abs(p2[1] - 0.0f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::pow", []()
              {
        using namespace Swarm;
        Point<float, 2> p1(2.0f);
        Point<float, 2> p2 = p1.pow(3);
        if (std::abs(p2[0] - 8.0f) > 1e-4f || std::abs(p2[1] - 8.0f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::norm1", []()
              {
        using namespace Swarm;
        Point<float, 3> p1([](size_t index)
                           { return static_cast<float>(index + 1); });
        float n1 = p1.norm1();
        if (std::abs(n1 - 6.0f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::squareNorm2", []()
              {
        using namespace Swarm;
        Point<float, 3> p1([](size_t index)
                           { return static_cast<float>(index + 1); });
        float n2 = p1.squareNorm2();
        if (std::abs(n2 - 14.0f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::norm2", []()
              {
        using namespace Swarm;
        Point<float, 3> p1([](size_t index)
                           { return static_cast<float>(index + 1); });
        float n2 = p1.norm2();
        if (std::abs(n2 - std::sqrt(14.0f)) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::abs", []()
              {
        using namespace Swarm;
        Point<float, 3> p1([](size_t index)
                           { float f = -1 * static_cast<float>((index + 1));
                            return f; });
        Point<float, 3> p2 = p1.abs();

        if (std::abs(p2[0] - 1.0f) > 1e-4f || std::abs(p2[1] - 2.0f) > 1e-4f || std::abs(p2[2] - 3.0f) > 1e-4f)
            return ERR;
        return PASS; });

    m.addTest("Test Point::isInsideBox", []()
              {
        using namespace Swarm;

        Point<float, 2> p1(1.0f);
        Point<float, 2> min(-2.0f);
        Point<float, 2> max(2.0f);
        
        if (!p1.isInsideBox(min, max))
            return ERR;

        Point<float, 2> p2(3.0f);
        if (p2.isInsideBox(min, max))
            return ERR;
        return PASS; });

    if (m.runAll() == ERR)
        return 1;

    return 0;
}