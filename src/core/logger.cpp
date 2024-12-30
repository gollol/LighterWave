#include <lightwave/logger.hpp>

namespace lightwave {
#ifdef DEBUG_PIXEL_LOGGER
thread_local t_debugPixelState debugPixelState;
#endif

Logger logger;
}
