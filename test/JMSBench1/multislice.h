// file: 'multislice.h'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file declares functions related to multislice calculations
// calling routines from the JMultiSlice library.
//
/* -----------------------------------------------------------------------

	This file is part of JMSBench1.

	JMSBench1 is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	JMSBench1 is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with JMSBench1.  If not, see <http://www.gnu.org/licenses/>.

----------------------------------------------------------------------- */

#pragma once

#include <thread>
#include "prm_main.h"
#include "JMultiSlice.h"

#ifndef __MS_WORKER_STATES
#define __MS_WORKER_STATES
#define MS_WS_IDLE			0
#define MS_WS_INIT			1
#define MS_WS_CALCULATE		2
#define MS_WS_OFFSET_PROBE	3
#define MS_WS_MULTISLICE	4
#define MS_WS_READOUT		5
#define MS_WS_ERROR			99
#define MS_WS_FINISHED		100
#define MS_WS_CANCELLED		98
#endif

class jmsworker {
public:

	// standard constructor
	jmsworker();

	// standard destructor
	~jmsworker();

	// member variables
public:

	bool b_wave_offset; // flags application of wave function shifts
	bool b_foc_avg; // apply explicit focal averaging
	bool b_active; // flags that the thread is still running
	bool b_cancel; // set to true if you want to cancel the calculation
	int n_state; // thread state flag
	int n_thread; // thread number of the application (-1: for no cpu)
	int n_gpu; // gpu id (-1: for no gpu)
	int whichcode; // flag for CPU or GPU code (_JMS_CODE_CPU or _JMS_CODE_GPU)
	int n_repeat; // number of repeats for accumulation
	int n_err; // error code
	unsigned int n_scan_idx; // scan pixel index
	unsigned int n_scan_ix; // scan pixel columns index
	unsigned int n_scan_iy; // scan pixel row index
	unsigned int n_acc_data; // data accumulation mode: 0 -> init with zero, 1 -> keep previous result
	unsigned int n_det_flags; // detection flags (see _JMS_DETECT_* in JMultislice.h)
	unsigned int n_fkern_num; // size of the focal convolution kernel (default: 7)
	float f_dx; // wave function shift along x in nm
	float f_dy; // wave function shift along y in nm
	float f_dz; // wave function shift along z in nm
	float f_fs; // focus spread in nm
	float f_fkw; // focal kernel relative width (default: 2.0f)
	float f_wgt; // weight of the result for accumulation
	float *pf_res_int; // result buffer for integrating detectors
	float *pf_res_img; // result buffer for image detectors
	float *pf_res_dif; // result buffer for diffraction detectors
	std::thread::id id_thread; // system cpu thread id
	std::string str_err; // error message
	CJMultiSlice *pjms; // pointer to multislice module
};


class jms_calc_queue
{
public:

	// standard constructor
	jms_calc_queue();

	// destructor
	~jms_calc_queue();

	// member variables

protected:

	unsigned int m_state; // queue state (0 = not ready, 1 = ready)
	jmsworker* m_q; // task queue
	size_t m_len_q; // number of queue tasks
	size_t m_num_tasks_open; // number of open tasks (those with state = 0)
	size_t m_num_tasks_solved; // number of open tasks (those with state = 3)
	size_t m_num_tasks_failed; // number of failed tasks (those with state = 5)
	mutable std::mutex guard;

	// member functions

public:

	// initializes the queue with a length of len_q tasks
	int init(size_t len_q);

	// sets the task data with index idx using data from jmsworker w
	int set_task(size_t idx, jmsworker task);

	// returns true if there are no open tasks to distribute
	bool empty(void);

	// return sthe number of open tasks
	size_t open(void);

	// returns the number of solved tasks
	size_t solved(void);

	// returns the number of failed tasks
	size_t failed(void);

	// returns the total number of tasks
	size_t total(void);

	// returns a pointer to a task for a requesting thread
	size_t request(void);

	// get specific task parameters
	int get_task(size_t idx, jmsworker& task_copy);

	// set task calculation state
	int set_task_calculate(size_t idx);

	// set task solved state
	int set_task_solved(size_t idx);

	// set task cancel state
	int set_task_cancel(size_t idx);

	// set task error state
	int set_task_error(size_t idx);
};


// prepare object transmission functions of a CJMultiSlice object
// assumes that phase gratings are ready in the sample.v_slc list of pprm
unsigned int __cdecl prepare_slices(prm_main *pprm, CJMultiSlice *pjms);

// prepare detector functions of a CJMultiSlice object
unsigned int __cdecl prepare_detectors(prm_main *pprm, CJMultiSlice *pjms);

// prepare calculation cores of a CJMultiSlice object
unsigned int __cdecl prepare_cores(prm_main *pprm, CJMultiSlice *pjms);

// prepare incident wave function of a CJMultiSlice object
// assumes that basic parameters of pjms are set: energy, grid, sampling
// also assumes that the calculation cores are initialized (memory allocations)
unsigned int __cdecl prepare_probe(prm_main *pprm, CJMultiSlice *pjms);


// Multislice calling interfaces

// single multislice call with one set of worker data
unsigned int __cdecl run_multislice(jmsworker* pw);

// cycling multislice calls with queued worker data
unsigned int __cdecl worker_multislice(int gpu_id, int cpu_id, prm_main* pprm, CJMultiSlice* pjms, jms_calc_queue* pq);


// runs a single thread STEM simulation on CPU or GPU
unsigned int __cdecl singlethread_stem(prm_main *pprm);

// runs a multi-thread STEM simulation on CPU and GPU
unsigned int __cdecl multithread_stem(prm_main *pprm);