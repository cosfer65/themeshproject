/* high resolution timer */
#ifndef __timer_h__
#define __timer_h__

namespace base_opengl {

    /**
     * @class gl_timer
     * @brief High-resolution timer for measuring elapsed time and system time.
     *
     * This class provides functionality to start, stop, reset, and advance a timer,
     * as well as retrieve absolute and elapsed times with high precision.
     * Platform-specific implementations are provided for Windows and POSIX systems.
     */
    class gl_timer
    {
    public:
        /**
         * @brief Constructs a new gl_timer object and initializes the timer.
         */
        gl_timer();

        /**
         * @brief Resets the timer to its initial state.
         */
        void reset();

        /**
         * @brief Starts or resumes the timer.
         */
        void start();

        /**
         * @brief Stops or pauses the timer.
         */
        void stop();

        /**
         * @brief Advances the timer by 0.1 seconds.
         */
        void advance();

        /**
         * @brief Gets the absolute system time in seconds.
         * @return The absolute system time as a double.
         */
        double get_absolute_time();

        /**
         * @brief Gets the current timer value in seconds.
         * @return The current timer value as a double.
         */
        double get_time();

        /**
         * @brief Gets the time elapsed since the last call to get_elapsed_time().
         * @return The elapsed time in seconds as a double.
         */
        double get_elapsed_time();

        /**
         * @brief Checks if the timer is currently stopped.
         * @return True if the timer is stopped, false otherwise.
         */
        bool is_stopped();

    protected:
        /// Indicates whether the timer is currently stopped.
        bool m_bTimerStopped;

        /// Indicates if QueryPerformanceFrequency is used.
        bool m_bUsingQPF;
        /// Ticks per second for QueryPerformanceFrequency.
        LONGLONG m_llQPFTicksPerSec;

        /// Time when the timer was stopped.
        LONGLONG m_llStopTime;
        /// Last elapsed time value.
        LONGLONG m_llLastElapsedTime;
        /// Base time for the timer.
        LONGLONG m_llBaseTime;
    };

    /**
     * @brief Gets a pointer to the global timer instance.
     * @return Pointer to the global gl_timer object.
     */
    gl_timer* get_global_timer();
}

#endif // __timer_h__
