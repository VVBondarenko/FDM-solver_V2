#include "gnuplointerface.h"

GnuploInterface::GnuploInterface()
{
    _pipe = nullptr;
}

GnuploInterface::~GnuploInterface()
{
    close();
}

bool GnuploInterface::isOpened() const
{
    return _pipe != nullptr;
}

void GnuploInterface::open()
{
    close();
    _pipe = popen("gnuplot", "w");
}

void GnuploInterface::flush()
{
    if (isOpened())
    {
        fflush(_pipe);
    }
}

void GnuploInterface::close()
{
    if (isOpened())
    {
        write("exit");
        flush();

        pclose(_pipe);
        _pipe = nullptr;
    }
}

void GnuploInterface::write(const char *line)
{
    if (isOpened() && line != nullptr && line[0] != '\0')
    {
        fprintf(_pipe, "%s\n", line);
    }
}

void GnuploInterface::write(const std::__cxx11::string &line)
{
    if (!line.empty())
    {
        write(line.c_str());
    }
}

void GnuploInterface::execute(const std::vector<std::__cxx11::string> &script)
{
    if (isOpened())
    {
        for (size_t i = 0; i < script.size(); i++)
        {
            write(script[i]);
            flush();
        }
    }
}
