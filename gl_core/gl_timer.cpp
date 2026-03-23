#define WIN32_LEAN_AND_MEAN
#define _CRT_SECURE_NO_WARNINGS
#define NOMINMAX
#include <windows.h>
#include "gl_graphics.h"
#include "gl_timer.h"

namespace base_opengl {

    /**
     * @brief Returns a pointer to the global timer instance.
     * 
     * This function provides access to a static instance of the gl_timer class,
     * ensuring a single timer is used throughout the application.
     * 
     * @return gl_timer* Pointer to the global timer instance.
     */
    gl_timer* get_global_timer()
    {
        static gl_timer timer;
        return &timer;
    }

    /**
     * @brief Constructs a gl_timer object and initializes timer state.
     * 
     * Initializes internal variables and attempts to use the high-resolution
     * performance counter if available.
     */
    gl_timer::gl_timer()
    {
        m_bTimerStopped = true;
        m_bUsingQPF = false;
        m_llQPFTicksPerSec = 0;
        m_llStopTime = 0;
        m_llLastElapsedTime = 0;
        m_llBaseTime = 0;

        LARGE_INTEGER qwTicksPerSec;
        m_bUsingQPF = (bool)(QueryPerformanceFrequency(&qwTicksPerSec) != 0);
        m_llQPFTicksPerSec = qwTicksPerSec.QuadPart;
    }

    /**
     * @brief Resets the timer to the current time.
     * 
     * Sets the base time and last elapsed time to the current value of the
     * performance counter. If the timer is stopped, uses the stop time.
     */
    void gl_timer::reset()
    {
        if (!m_bUsingQPF)
            return;

        LARGE_INTEGER qwTime;
        if (m_llStopTime != 0)
            qwTime.QuadPart = m_llStopTime;
        else
            QueryPerformanceCounter(&qwTime);

        m_llBaseTime = qwTime.QuadPart;
        m_llLastElapsedTime = qwTime.QuadPart;
        m_llStopTime = 0;
        m_bTimerStopped = FALSE;
    }

    /**
     * @brief Starts or resumes the timer.
     * 
     * If the timer was previously stopped, adjusts the base time to account for
     * the paused duration and resumes timing.
     */
    void gl_timer::start()
    {
        if (!m_bUsingQPF)
            return;

        LARGE_INTEGER qwTime;
        QueryPerformanceCounter(&qwTime);

        if (m_bTimerStopped)
            m_llBaseTime += qwTime.QuadPart - m_llStopTime;
        m_llStopTime = 0;
        m_llLastElapsedTime = qwTime.QuadPart;
        m_bTimerStopped = FALSE;
    }

    /**
     * @brief Stops or pauses the timer.
     * 
     * Records the current time as the stop time and marks the timer as stopped.
     */
    void gl_timer::stop()
    {
        if (!m_bUsingQPF)
            return;

        if (!m_bTimerStopped)
        {
            LARGE_INTEGER qwTime;
            if (m_llStopTime != 0)
                qwTime.QuadPart = m_llStopTime;
            else
                QueryPerformanceCounter(&qwTime);

            m_llStopTime = qwTime.QuadPart;
            m_llLastElapsedTime = qwTime.QuadPart;
            m_bTimerStopped = TRUE;
        }
    }

    /**
     * @brief Advances the timer by a fixed interval (0.1 seconds).
     * 
     * Useful for simulating time progression when the timer is stopped.
     */
    void gl_timer::advance()
    {
        if (!m_bUsingQPF)
            return;
        m_llStopTime += m_llQPFTicksPerSec / 10;
    }

    /**
     * @brief Gets the absolute system time in seconds.
     * 
     * Returns the current value of the performance counter in seconds.
     * 
     * @return double Absolute time in seconds, or -1.0 if unavailable.
     */
    double gl_timer::get_absolute_time()
    {
        if (!m_bUsingQPF)
            return -1.0;

        LARGE_INTEGER qwTime;
        if (m_llStopTime != 0)
            qwTime.QuadPart = m_llStopTime;
        else
            QueryPerformanceCounter(&qwTime);

        double fTime = qwTime.QuadPart / (double)m_llQPFTicksPerSec;

        return fTime;
    }

    /**
     * @brief Gets the application time in seconds since the last reset.
     * 
     * Returns the elapsed time since the timer was last reset or started.
     * 
     * @return double Application time in seconds, or -1.0 if unavailable.
     */
    double gl_timer::get_time()
    {
        if (!m_bUsingQPF)
            return -1.0;

        LARGE_INTEGER qwTime;
        if (m_llStopTime != 0)
            qwTime.QuadPart = m_llStopTime;
        else
            QueryPerformanceCounter(&qwTime);

        double fAppTime = (double)(qwTime.QuadPart - m_llBaseTime) / (double)m_llQPFTicksPerSec;

        return fAppTime;
    }

    /**
     * @brief Gets the time elapsed since the last call to get_elapsed_time().
     * 
     * Updates the last elapsed time to the current time.
     * 
     * @return double Elapsed time in seconds, or -1.0 if unavailable.
     */
    double gl_timer::get_elapsed_time()
    {
        if (!m_bUsingQPF)
            return -1.0;

        LARGE_INTEGER qwTime;
        if (m_llStopTime != 0)
            qwTime.QuadPart = m_llStopTime;
        else
            QueryPerformanceCounter(&qwTime);

        double fElapsedTime = (double)(qwTime.QuadPart - m_llLastElapsedTime) / (double)m_llQPFTicksPerSec;
        m_llLastElapsedTime = qwTime.QuadPart;

        return fElapsedTime;
    }

    /**
     * @brief Checks if the timer is currently stopped.
     * 
     * @return true if the timer is stopped, false otherwise.
     */
    bool gl_timer::is_stopped()
    {
        return m_bTimerStopped;
    }
}