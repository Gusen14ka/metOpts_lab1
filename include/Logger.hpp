#pragma once

#include <fstream>
#include <string>

class Logger{
public:
    enum class Level {
        INFO,
        WARNING,
        ERR
    };

    static Logger& instance();
    ~Logger();

    void info(std::string const & message, std::string const & component);
    void warning(std::string const & message, std::string const & component);
    void error(std::string const & message, std::string const & component);

    void set_log_to_console(bool enable);
    void set_log_to_file(bool enable);
    void set_log_file(std::string const & path);

    Logger(Logger const &) = delete;
    Logger& operator=(Logger const &) = delete;
    Logger(Logger const &&) = delete;
    Logger& operator=(Logger const &&) = delete;

private:
    Logger() = default;
    std::string level_to_str(Level lv);
    void log(Level lv, std::string const & message, std::string const & component);

    std::string logFilePath_;
    std::ofstream logFile_;
    bool consoleLog_ = true;
    bool fileLog_ = false;
};