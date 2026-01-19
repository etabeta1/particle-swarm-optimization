#ifndef _TESTS_HPP_
#define _TESTS_HPP_

#include <string>
#include <functional>
#include <vector>
#include <iostream>

enum Status
{
    PASS,
    ERR
};

struct Test
{
    std::string name;
    std::function<Status()> body;
};

class TestManager
{
private:
    std::vector<Test> tests;

public:
    TestManager() {}

    void addTest(const std::string &name, const std::function<Status()> &body)
    {
        Test t{name, body};
        tests.emplace_back(t);
    }

    Status runAll()
    {
        size_t pass = 0, err = 0;

        std::cout << "Running " << tests.size() << " test(s)" << std::endl;

        for (Test &test : tests)
        {
            std::cout << "Running test \"" << test.name << "\"... ";
            Status r = test.body();

            switch (r)
            {
            case PASS:
                std::cout << "[ \033[32mPASS\033[0m ]";
                pass++;
                break;

            default:
                std::cout << "[ \033[31mERR\033[0m ]";
                err++;
                break;
            }

            std::cout << std::endl;
        }

        std::cout << "Final results:\n  PASS:\t" << pass << "\n  ERR:\t" << err << std::endl;

        return err == 0 ? PASS : ERR;
    }
};

#endif