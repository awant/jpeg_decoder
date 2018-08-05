#pragma once

#include <iostream>

#ifdef LOGGING_ENABLED
#define DEBUG(msg) std::cout << "debug: " << msg << " @" << __FILE__ << ":" << __LINE__ << std::endl;
#define INFO(msg) std::cout << "info: " << msg << " @" << __FILE__ << ":" << __LINE__ << std::endl;
#define WARNING(msg) std::cout << "warning: " << msg << " @" << __FILE__ << ":" << __LINE__ << std::endl;
#define ERROR(msg) std::cerr << "error: " << msg << " @" << __FILE__ << ":" << __LINE__ << std::endl;
#else
#define LOG(level, msg)
#endif
