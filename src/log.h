#ifndef YICPPLIB_LOG_H
#define YICPPLIB_LOG_H

#include <iostream>

namespace YiCppLib {
    namespace Log {
        typedef struct {
            unsigned int loglvl; 
            const char * logname;
        } LogLvlT;
        static std::ostream bitBucket(0);
    }
}

#define LOG(loglvl) (std::cerr<<"["<<loglvl.logname<<"] ")

#define LOGGER(maxLogLvl) ([](YiCppLib::Log::LogLvlT loglvl) -> std::ostream& { \
        return loglvl.loglvl >= maxLogLvl.loglvl ? LOG(loglvl) : YiCppLib::Log::bitBucket; })

#define LOGLV_DEBUG ((YiCppLib::Log::LogLvlT){.loglvl=0, .logname="DEBUG"})
#define LOGLV_INFO ((YiCppLib::Log::LogLvlT){.loglvl=1, .logname="INFO"})
#define LOGLV_WARN ((YiCppLib::Log::LogLvlT){.loglvl=2, .logname="WARN"})
#define LOGLV_ERR ((YiCppLib::Log::LogLvlT){.loglvl=3, .logname="ERROR"})

#ifndef LOGLVL
auto logger = LOGGER(LOGLV_WARN);
#else
auto logger = LOGGER(LOGLVL);
#endif

#endif
