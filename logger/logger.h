#pragma once

#include <iostream>
#include <string>
#include <ctime>
#include <cstdlib>

namespace da {
    struct nullstream: std::ostream {
        nullstream(): std::ostream(nullptr) {}
    };

    template <typename T>
    nullstream &operator<<(nullstream &o, T const & x) { return o;}

    nullstream  __nullstream;

    class LogMessage {
    public:
        explicit LogMessage(const std::string& l)
                : level_(l), ofs_(enable_ ? (l == "ERROR" ? std::cerr : std::cout) : __nullstream) {
            stream() << "[" << level_ << "]\t";
        }

        explicit LogMessage(std::ostream &o)
                : level_("ERROR"), ofs_(o) {
            stream() << "[" << level_ << "]\t";
        }

        inline std::ostream &stream() {
            return ofs_;
        }

        ~LogMessage() {
            stream() << std::endl;
        }

        static void Enable(bool enable) {
            enable_ = enable;
        }

    private:
        std::string level_;
        std::ostream &ofs_;
        static bool enable_;
    };
    bool LogMessage::enable_ = true;
}

#define LOG(type)        da::LogMessage(#type).stream()
#define LOG_ERROR        LOG(ERROR)
#define LOG_DEBUG        LOG(DEBUG)
