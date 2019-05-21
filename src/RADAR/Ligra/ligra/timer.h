// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef TIMER_H_
#define TIMER_H_

#include <sys/time.h>
#include <ctime>
#include <chrono>


/*
GAP Benchmark Suite
Class:  Timer
Author: Scott Beamer

Simple timer that wraps gettimeofday
*/

#define COMPILER_BARRIER() asm volatile("" ::: "memory")
 
class Timer {
 public:
  Timer() { m_elapsedTime = 0.0;}

  void Start() {
    COMPILER_BARRIER(); 
    start_time_ = std::chrono::steady_clock::now();
    COMPILER_BARRIER(); 
  }

  void Stop() {
    COMPILER_BARRIER(); 
    stop_time_ = std::chrono::steady_clock::now();
    COMPILER_BARRIER(); 
    auto diff = stop_time_ - start_time_;
    m_elapsedTime += std::chrono::duration<double>(diff).count(); 
  }

  double Seconds() const {
    return m_elapsedTime; 
  }

  double Millisecs() const {
    return m_elapsedTime * 1000; 
  }

  double Microsecs() const {
    return m_elapsedTime * 1000000; 
  }
  
  void PrintTime(const std::string &s, double seconds) {
    printf("%-21s%3lf\n", (s + ":").c_str(), seconds);
  }

 private:
  std::chrono::steady_clock::time_point start_time_;
  std::chrono::steady_clock::time_point stop_time_;
  double m_elapsedTime;
};

// Times op's execution using the timer t
#define TIME_OP(t, op) { t.Start(); (op); t.Stop(); }

#endif  // TIMER_H_
