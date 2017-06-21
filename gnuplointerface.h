#ifndef GNUPLOINTERFACE_H
#define GNUPLOINTERFACE_H

/*
 * GNU Plot Interface by Can Erhan
 * copy-pasted from
 *  https://github.com/ccerhan/gnuplot-cpp-interface
 *
 * Changes:
 *
 */

#include <cstdio>
#include <cstdlib>
#include <ostream>
#include <string>
#include <vector>

class GnuploInterface
{
private:
    FILE* _pipe;

public:
    GnuploInterface();
    virtual ~GnuploInterface();

    bool isOpened() const;

    void open();
    void flush();
    void close();
    void write(const char *line);
    void write(const std::string &line);
    void execute(const std::vector<std::string> &script);
};

#endif // GNUPLOINTERFACE_H
