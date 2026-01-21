#ifndef _TESTS_HPP_
#define _TESTS_HPP_

#include <string>
#include <functional>
#include <vector>
#include <iostream>

/**
 * @brief Status of a test
 */
enum Status
{
    /**
     * @brief Test passed
     */
    PASS,
    /**
     * @brief Test failed
     */
    ERR
};

/**
 * @brief Structure representing a test
 */
struct Test
{
    /**
     * @brief Name of the test
     */
    std::string name;

    /**
     * @brief Body of the test
     */
    std::function<Status()> body;
};

/**
 * @brief Manager for tests
 */
class TestManager
{
private:
    /**
     * @brief List of tests
     */
    std::vector<Test> tests;

public:
    /**
     * @brief Construct a new Test Manager object
     */
    TestManager() {}

    /**
     * @brief Add a test to the manager
     * @param name Name of the test
     * @param body Body of the test
     */
    void addTest(const std::string &name, const std::function<Status()> &body)
    {
        Test t{name, body};
        tests.emplace_back(t);
    }

    /**
     * @brief Run all tests
     * @return Status Overall status (`PASS` if all tests passed, `ERR` otherwise)
     */
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