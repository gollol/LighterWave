#include <lightwave/camera.hpp>
#include <lightwave/integrator.hpp>
#include <lightwave/parallel.hpp>

#include <algorithm>
#include <chrono>

#include <lightwave/iterators.hpp>
#include <lightwave/streaming.hpp>

namespace lightwave {

void SamplingIntegrator::execute() {
    if (!m_image) {
        lightwave_throw(
            "<integrator /> needs an <image /> child to render into!");
    }

    const Vector2i resolution = m_scene->camera()->resolution();
    m_image->initialize(resolution);

    const float norm = 1.0f / m_sampler->samplesPerPixel();

    Streaming stream{ *m_image };
    ProgressReporter progress{ resolution.product() };
    for_each_parallel(BlockSpiral(resolution, Vector2i(64)), [&](auto block) {
        auto sampler = m_sampler->clone();

        #ifdef DEBUG_PIXEL_LOGGER
        debugPixelState.active = true;
        #endif

        for (auto pixel : block) {
            #ifdef DEBUG_PIXEL_LOGGER
            debugPixelState.x = pixel.x();
            debugPixelState.y = pixel.y();
            #endif

            Color sum;  
            for (int sample = 0; sample < m_sampler->samplesPerPixel();
                 sample++) {
                #ifdef DEBUG_PIXEL_LOGGER
                debugPixelState.sample = sample;
                #endif

                Context cont = {.rng = &*sampler, .scene = &*m_scene};
                sampler->seed(pixel, sample);
                auto cameraSample = m_scene->camera()->sample(pixel, *sampler);
                sum += cameraSample.weight * Li(cameraSample.ray, *sampler, cont);
            }
            m_image->get(pixel) = norm * sum;
        }

        #ifdef DEBUG_PIXEL_LOGGER
        debugPixelState.active = false;
        #endif

        progress += block.diagonal().product();
        stream.updateBlock(block);
    });
    progress.finish();

    m_image->save();
}

} // namespace lightwave