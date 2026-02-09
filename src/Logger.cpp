#include "Logger.hpp"
#include <iostream>
#include <string>

Logger& Logger::instance(){
    static Logger logger;
    return logger;
}

Logger::~Logger(){
    if(logFile_.is_open()){
        logFile_.close();
    }
}

void Logger::info(std::string const & message, std::string const & component){
    log(Level::INFO, message, component);
}

void Logger::warning(std::string const & message, std::string const & component){
    log(Level::WARNING, message, component);
}

void Logger::error(std::string const & message, std::string const & component){
    log(Level::ERR, message, component);
}

void Logger::set_log_to_console(bool enable){
    consoleLog_ = enable;
}

void Logger::set_log_to_file(bool enable){
    fileLog_ = enable;
    if(enable && !logFile_.is_open() && !logFilePath_.empty()){
        logFile_.open(logFilePath_);
        if(!logFile_.is_open() || logFile_.fail()){
            fileLog_ = false;
            logFile_.close();
        }
    }
}

void Logger::set_log_file(std::string const & path){
    logFilePath_ = path;
}

std::string Logger::level_to_str(Level lv){
    switch (lv) {
        case Level::INFO: return "INFO";
        case Level::WARNING: return "WARNING";
        case Level::ERR: return "ERROR";
        default: return "UNKNOWN";
    }
}

void Logger::log(Level lv, std::string const & message, std::string const & component){
    std::string output;
    output = "[" + level_to_str(lv) + "] ";
    if(!component.empty()) output +=  "[" + component + "] ";
    output += message;

    if(consoleLog_){
        if(lv == Level::ERR){
            std::cerr << output << std::endl;
        }
        else{
            std::cout << output << std::endl;
        }
    }
    if(fileLog_ && logFile_.is_open()){
        logFile_ << output << std::endl;
    }
}