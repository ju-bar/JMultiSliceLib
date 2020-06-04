// file: 'multislice.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains functions related to multislice calculations
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

#include "multislice.h"
#include "string_format.h"


jmsworker::jmsworker()
{
	b_active = false;
	b_cancel = false;
	b_foc_avg = false;
	b_wave_offset = true;
	whichcode = 0;
	f_dx = 0.f;
	f_dy = 0.f;
	f_dz = 0.f;
	f_fkw = 2.f;
	f_fs = 0.f;
	f_wgt = 0.f;
	id_thread = std::thread::id();
	n_acc_data = _JMS_ACCMODE_NONE;
	n_det_flags = _JMS_DETECT_NONE;
	n_err = 0;
	n_fkern_num = 7;
	n_gpu = -1;
	n_repeat = 1;
	n_scan_idx = 0;
	n_scan_ix = 0;
	n_scan_iy = 0;
	n_state = MS_WS_IDLE;
	n_thread = -1;
	pf_res_dif = NULL;
	pf_res_img = NULL;
	pf_res_int = NULL;
	pjms = NULL;
}

jmsworker::~jmsworker()
{

}

// ----------------------------------------------------------------------------
//
// class jms_calc_queue
//
// ----------------------------------------------------------------------------

jms_calc_queue::jms_calc_queue()
{
	m_state = 0;
	m_q = 0;
	m_len_q = 0;
	m_num_tasks_open = 0;
	m_num_tasks_solved = 0;
	m_num_tasks_failed = 0;
}

jms_calc_queue::~jms_calc_queue()
{
	if (NULL != m_q) { delete[] m_q; }
}

int jms_calc_queue::init(size_t len_q)
{
	int nerr = 0;

	if (NULL != m_q) {
		delete[] m_q;
		m_q = NULL;
		m_state = 0;
		m_len_q = 0;
		m_num_tasks_open = 0;
		m_num_tasks_solved = 0;
		m_num_tasks_failed = 0;
	}

	if (len_q == 0) {
		// finish without tasks
		goto _exit;
	}

	m_q = new jmsworker[len_q];
	if (NULL == m_q) {
		nerr = 1;
		std::cerr << "Error (jms_calc_queue::init): Failed to initialize queue of length " << len_q << "." << std::endl;
		goto _exit;
	}
	m_num_tasks_solved = 0;
	m_num_tasks_failed = 0;
	m_num_tasks_open = len_q;
	m_len_q = len_q;
	m_state = 1;

_exit:
	return nerr;
}

int jms_calc_queue::set_task(size_t idx, jmsworker task)
{
	int nerr = 0;
	if (idx >= m_len_q) {
		nerr = 1;
		std::cerr << "Error (jms_calc_queue::set_task): Invalid task index (" << idx << ") for current queue of length " << m_len_q << "." << std::endl;
		goto _exit;
	}
	m_q[idx] = task;
_exit:
	return nerr;
}

bool jms_calc_queue::empty(void)
{
	return (0 == open());
}

size_t jms_calc_queue::open(void)
{
	return m_num_tasks_open;
}

size_t jms_calc_queue::solved(void)
{
	return m_num_tasks_solved;
}

size_t jms_calc_queue::failed(void)
{
	return m_num_tasks_failed;
}

size_t jms_calc_queue::total(void)
{
	return m_len_q;
}

size_t jms_calc_queue::request(void)
{
	size_t itask = 0;
	if (!empty()) {
		std::lock_guard<std::mutex> lock(guard);
		for (itask = 0; itask < m_len_q; itask++) {
			if (m_q[itask].n_state == MS_WS_IDLE) {
				m_q[itask].n_state = MS_WS_INIT;
				m_num_tasks_open--;
				break;
			}
		}
		return itask;
	}
	return NULL;
}

// get specific task parameters
int jms_calc_queue::get_task(size_t idx, jmsworker& task_copy)
{
	if (idx < m_len_q) {
		if (m_q[idx].n_state == MS_WS_INIT) {
			task_copy = m_q[idx];
			return 0;
		}
	}
	return 1;
}

// set task calculation state
int jms_calc_queue::set_task_calculate(size_t idx)
{
	if (idx < m_len_q) {
		if (m_q[idx].n_state == MS_WS_INIT) {
			m_q[idx].n_state = MS_WS_CALCULATE;
			return 0;
		}
	}
	return 1;
}

// set task solved state
int jms_calc_queue::set_task_solved(size_t idx)
{
	if (idx < m_len_q) {
		if (m_q[idx].n_state == MS_WS_CALCULATE) {
			m_q[idx].n_state = MS_WS_FINISHED;
			m_num_tasks_solved++;
			return 0;
		}
	}
	return 1;
}

// set task cancel state
int jms_calc_queue::set_task_cancel(size_t idx)
{
	if (idx < m_len_q) {
		if (m_q[idx].n_state == MS_WS_CALCULATE) {
			m_q[idx].n_state = MS_WS_CANCELLED;
			return 0;
		}
	}
	return 1;
}

// set task error state
int jms_calc_queue::set_task_error(size_t idx)
{
	if (idx < m_len_q) {
		if (m_q[idx].n_state == MS_WS_CALCULATE) {
			m_q[idx].n_state = MS_WS_ERROR;
			m_num_tasks_failed++;
			return 0;
		}
	}
	return 1;
}







unsigned int __cdecl prepare_slices(prm_main *pprm, CJMultiSlice *pjms)
{
	unsigned int nerr = 0;
	int njmserr = 0;
	int whichcode = 0;
	int *proidx = NULL;

	if (NULL == pprm || NULL == pjms) {
		std::cerr << "Error: (prepare_slices) called with invalid arguments." << std::endl;
		nerr = 1;
		goto _exit;
	}

	if (pprm->cpu_num > 0) whichcode += (int)_JMS_CODE_CPU;
	if (pprm->gpu_id >= 0) whichcode += (int)_JMS_CODE_GPU;
	if (0 == whichcode) goto _exit;

	if (pprm->btalk) {
		std::cout << std::endl;
		std::cout << "  Preparing object transmission functions ..." << std::endl;
	}

	// setup phase grating memory in JMS
	njmserr = pprm->sample.setup_pgr_jms(pjms, whichcode);
	if (0 < njmserr) {
		std::cerr << "Error: (prepare_slices): failed to setup phase gratings in multislice module (" << njmserr << ")." << std::endl;
		nerr = 2;
		goto _exit;
	}
		
	// setup object slice sequence in JMS
	njmserr = pprm->sample.setup_object_slice_seq_jms(pjms);
	if (0 < njmserr) {
		std::cerr << "Error: (prepare_slices): failed to setup object slice sequence in multislice module (" << njmserr << ")." << std::endl;
		nerr = 4;
		goto _exit;
	}

	// prepare propagator data
	njmserr = pprm->sample.prepare_pro(pjms, whichcode);
	if (0 < njmserr) {
		std::cerr << "Error: (prepare_slices): failed to calculate propagators by multislice module (" << njmserr << ")." << std::endl;
		nerr = 5;
		goto _exit;
	}

_exit:
	return nerr;
}

unsigned int __cdecl prepare_detectors(prm_main *pprm, CJMultiSlice *pjms)
{
	unsigned int nerr = 0;
	int num_det = 0; // number of STEM detectors
	int det_flags = 0; // detection flags
	int njmserr = 0;
	int whichcode = 0;
	int idet = 0;

	if (NULL == pprm || NULL == pjms) {
		std::cerr << "Error: (prepare_detectors) called with invalid arguments." << std::endl;
		nerr = 1;
		goto _exit;
	}

	if (pprm->cpu_num > 0) whichcode += (int)_JMS_CODE_CPU;
	if (pprm->gpu_id >= 0) whichcode += (int)_JMS_CODE_GPU;
	if (0 == whichcode) goto _exit;

	if (pprm->btalk) {
		std::cout << std::endl;
		std::cout << "  Preparing detector functions ..." << std::endl;
	}

	num_det = (int)pprm->detector.v_annular.size();
	det_flags = pprm->detector.get_jms_flags();
	njmserr = pjms->DetectorSetup(whichcode, pprm->sample.v_slc_det, num_det, det_flags, pprm->cpu_num);
	if (0 < njmserr) {
		std::cerr << "Error: (prepare_detectors): failed to prepare detector arrays (" << njmserr << ")." << std::endl;
		nerr = 2;
		goto _exit;
	}

	if (pprm->detector.b_annular && (0 < num_det)) {
		// prepare detector functions of all integrating detectors
		for (idet = 0; idet < num_det; idet++) {
			njmserr = pprm->detector.v_annular[idet].calculate_func_jms(pjms);
			if (0 < njmserr) {
				std::cerr << "Error: (prepare_detectors): failed to calculate detector #" << idet << " arrays (" << njmserr << ")." << std::endl;
				nerr = 100 + idet;
				goto _exit;
			}
			njmserr = pprm->detector.v_annular[idet].set_func_jms(pjms, whichcode, idet);
			if (0 < njmserr) {
				std::cerr << "Error: (prepare_detectors): failed to transfer detector function #" << idet << " to multislice module (" << njmserr << ")." << std::endl;
				nerr = 200 + idet;
				goto _exit;
			}
		}
	}

_exit:
	return nerr;
}

unsigned int __cdecl prepare_cores(prm_main *pprm, CJMultiSlice *pjms)
{
	unsigned int nerr = 0;
	int njmserr = 0;
	int whichcode = 0;

	if (NULL == pprm || NULL == pjms) {
		std::cerr << "Error: (prepare_cores) called with invalid arguments." << std::endl;
		nerr = 1;
		goto _exit;
	}

	if (pprm->cpu_num > 0) whichcode += (int)_JMS_CODE_CPU;
	if (pprm->gpu_id >= 0) whichcode += (int)_JMS_CODE_GPU;
	if (0 == whichcode) goto _exit;

	if (pprm->btalk) {
		std::cout << std::endl;
		std::cout << "  Initializing multislice cores ..." << std::endl;
		if (whichcode & (int)_JMS_CODE_GPU) {
			std::cout << "  - GPU core on device #" << pprm->gpu_id << "." << std::endl;
		}
		if (whichcode & (int)_JMS_CODE_CPU) {
			std::cout << "  - " << pprm->cpu_num << " CPU cores." << std::endl;
		}
	}

	njmserr = pjms->InitCore(whichcode, pprm->cpu_num);
	if (0 < njmserr) {
		nerr = 2;
		std::cerr << "Error: (prepare_cores) failed to initialize the multislice core (" << njmserr << ")." << std::endl;
		goto _exit;
	}

_exit:
	return nerr;
}

unsigned int __cdecl prepare_probe(prm_main *pprm, CJMultiSlice *pjms)
{
	unsigned int nerr = 0, ierr = 0;
	int njmserr = 0;
	fcmplx *wave = 0; // wave function buffer
	size_t num_pix = 0; // number of wave function pixels
	size_t sz_wave = 0; // size of wave function buffer in bytes
	int whichcode = 0;
	int nx = 0, ny = 0;
	CJProbeParams jpp;

	if (NULL == pprm || NULL == pjms) {
		std::cerr << "Error: (prepare_probe) called with invalid arguments." << std::endl;
		nerr = 1;
		goto _exit;
	}

	if (pprm->cpu_num > 0) whichcode += (int)_JMS_CODE_CPU;
	if (pprm->gpu_id >= 0) whichcode += (int)_JMS_CODE_GPU;
	if (0 == whichcode) goto _exit;

	if (pprm->btalk) {
		std::cout << std::endl;
		std::cout << "  Calculating incident probe wave function ..." << std::endl;
	}

	pjms->GetGridSize(nx, ny);
	num_pix = (size_t)nx * ny;
	sz_wave = sizeof(fcmplx) * num_pix;

	if (0 == num_pix) {
		std::cerr << "Error: (prepare_probe) zero grid sampling not supported." << std::endl;
		nerr = 2;
		goto _exit;
	}

	wave = (fcmplx*)malloc(sz_wave);
	if (NULL == wave) {
		std::cerr << "Error: (prepare_probe) failed to allocate memory for probe wave function." << std::endl;
		nerr = 3;
		goto _exit;
	}
	
	ierr = pprm->probe.get_jpg_params(&jpp);
	if (0 < ierr) {
		std::cerr << "Error: (prepare_probe) failed to allocate memory for probe wave function." << std::endl;
		nerr = 4;
		goto _exit;
	}
	
	// insert switch STEM / TEM here
	njmserr = pjms->CalculateProbeWaveFourier(&jpp, wave); // calculate STEM probe
	if (0 < njmserr) {
		std::cerr << "Error: (prepare_probe) failed to calculate probe wave function (" << njmserr << ")." << std::endl;
		nerr = 5;
		goto _exit;
	}

	njmserr = pjms->SetIncidentWave(whichcode, wave);
	if (0 < njmserr) {
		std::cerr << "Error: (prepare_probe) failed to transfer wavefunction to multislice module (" << njmserr << ")." << std::endl;
		nerr = 6;
		goto _exit;
	}

_exit:
	if (NULL != wave) {
		free(wave);
		wave = NULL;
	}
	return nerr;
}



unsigned int __cdecl run_multislice(jmsworker* pw)
{
	int nerr = 0; // routine error code
	int njmserr = 0; // jms error code
	int nrepeat = 1; // number of multislice repeats
	int nznum = 1; // number of focal steps
	int iz = 0; // focal index
	int ir = 0; // repeat index
	int ithread = 0; // thread index
	int nacc = _JMS_ACCMODE_NONE;
	float f_wgt = 1.f; // weight of a run;
	float f_wgt_cur = 1.f; // temp. weight
	float f_wgt_tot = 0.f; // total weight
	float fk_prm = 0.f; // focal kernel parameter
	float fk_0 = 0.f; // focal kernel offset
	float fk_step = 0.f; // focal kernel step
	float f_x = 0.f, f_y = 0.f, f_z = 0.f; // positions in nm
	std::string scode = "GPU";
	if (NULL == pw) {
		nerr = 1;
		goto _exit;
	}
	// worker state init
	pw->b_active = true;
	pw->b_cancel = false;
	pw->n_state = MS_WS_CALCULATE;
	pw->n_err = 0;
	pw->str_err = "";
	if (NULL == pw->pjms) {
		nerr = 2;
		pw->n_err = nerr;
		pw->str_err = "Invalid pointer to JMultiSlice module.";
		goto _exit;
	}
	if (pw->whichcode == _JMS_CODE_CPU) {
		scode = "CPU";
		ithread = pw->n_thread;
	}
	nacc = pw->n_acc_data;
	if (nacc == _JMS_ACCMODE_NONE) {
		pw->f_wgt = 0.f; // reset job weight
	}
	// focal averaging init
	if (pw->b_foc_avg) {
		nznum = pw->n_fkern_num;
		if (1 >= nznum) { // turn focal averaging off
			nznum = 1;
			pw->b_foc_avg = false;
		}
		if (fabsf(pw->f_fs) == 0.0f || fabsf(pw->f_fkw) == 0.0f) { // turn focal averaging off
			nznum = 1;
			pw->b_foc_avg = false;
		}
	}
	// repetition init
	nrepeat = pw->n_repeat;
	if (nrepeat < 1) nrepeat = 1;
	f_wgt = 1.f / (float)nrepeat; // init weight for each pass
	if (nznum > 1) { // make focal kernel size a factor of nrepeat
		nrepeat = nznum * (1 + (nrepeat - (nrepeat%nznum)) / nznum);
	}
	// finish focal kernel setup
	if (pw->b_foc_avg) { 
		fk_step = 2.0f * fabsf(pw->f_fkw) * fabsf(pw->f_fs) / (float)(nznum - 1); // focal setp size
		fk_0 = -1.f * fabsf(pw->f_fkw) * fabsf(pw->f_fs); // focal offset
		fk_prm = -1.0f / (fabsf(pw->f_fs) * fabsf(pw->f_fs)); // initialize exponent pre-factor 1/fs^2
		float fsum = 0.f;
		for (iz = 0; iz < nznum; iz++) {
			f_z = fk_0 + fk_step * (float)iz;
			fsum += exp(fk_prm * f_z * f_z);
		}
		f_wgt = (float)nznum / fsum / (float)nrepeat; // rescaled weight to single repeat step
	}

	//
	// running the multislice
	//
	f_wgt_tot = 0.f; // reset accumulation of weights
	iz = 0; // reset focal loop index
	f_x = pw->f_dx; // set probe x shift
	f_y = pw->f_dy; // set probe y shift
	//
	for (ir = 0; ir < nrepeat; ir++) {
		//
		if (pw->b_cancel) { // cancel check before new run
			goto _exit;
		}
		//
		f_wgt_cur = f_wgt; // preset the weight for the current repeat
		//
		// *** probe offset
		//
		if (pw->b_wave_offset) { // apply wave function offsets
			f_z = 0.f; // init probe z shift
			if (pw->b_foc_avg) { // add focal kernel z shift
				f_z = fk_0 + fk_step * (float)iz;
				f_wgt_cur *= exp(fk_prm * f_z * f_z); // update weight of current pass for focal kernel step
				iz++; // increment focal kernel step
				if (iz >= nznum) { // reset focal kernel step
					iz = 0;
				}
			}
			f_z += pw->f_dz; // add mean defocus
			njmserr = pw->pjms->OffsetIncomingWave(pw->whichcode, f_x, f_y, f_z, ithread); // call jms to offset the wf
			if (0 < njmserr) { // handle error
				nerr = 3;
				pw->n_err = nerr;
				pw->str_err = "Wave function offset failed on " + scode + " (" + format("%d", njmserr) + ")";
				goto _exit;
			}
			//
			if (pw->b_cancel) { // cancel check after probe offset
				goto _exit;
			}
		}
		//
		// *** multislice
		//
		if (pw->whichcode == _JMS_CODE_CPU) {
			njmserr = pw->pjms->CPUMultislice(0, nacc, f_wgt_cur, ithread);
		}
		else {
			njmserr = pw->pjms->GPUMultislice(0, nacc, f_wgt_cur);
		}
		f_wgt_tot += f_wgt_cur;
		pw->f_wgt += f_wgt_cur;
		if (0 < njmserr) { // handle error
			nerr = 4;
			pw->n_err = nerr;
			pw->str_err = "Multislice failed on " + scode + " (" + format("%d", njmserr) + ")";
			goto _exit;
		}
		nacc = _JMS_ACCMODE_INTEGRATE; // set local flag to integrating mode after the first pass
	}
	//
	// *** get the results
	//
	// The following code assumes, that the results collected in JMS are already scaled properly.
	// The weights used above should in total give 1.0. See f_wgt_tot ! At this point it should read 1.0
	//
	if (pw->b_cancel) { // cancel check before result acquisition
		goto _exit;
	}
	//
	if (pw->pjms->GetDetetionSliceNum() > 0) { // there is detection in general
		if ((pw->n_det_flags & (unsigned int)_JMS_DETECT_INTEGRATED) && (pw->pjms->GetDetNum() > 0)) {
			// integerating detector readout
			njmserr = pw->pjms->GetResult(pw->whichcode, _JMS_DETECT_INTEGRATED, pw->pf_res_int, ithread);
			if (0 < njmserr) {
				nerr = 11;
				pw->n_err = nerr;
				pw->str_err = "Retrieving integrating detector data failed on " + scode + " (" + format("%d", njmserr) + ")";
				goto _exit;
			}
		}
		if (pw->n_det_flags & (unsigned int)_JMS_DETECT_IMAGE) {
			// image detector readout
			njmserr = pw->pjms->GetResult(pw->whichcode, _JMS_DETECT_IMAGE, pw->pf_res_img, ithread);
			if (0 < njmserr) {
				nerr = 12;
				pw->n_err = nerr;
				pw->str_err = "Retrieving image data failed on " + scode + " (" + format("%d", njmserr) + ")";
				goto _exit;
			}
		}
		if (pw->n_det_flags & (unsigned int)_JMS_DETECT_DIFFRACTION) {
			// image detector readout
			njmserr = pw->pjms->GetResult(pw->whichcode, _JMS_DETECT_DIFFRACTION, pw->pf_res_dif, ithread);
			if (0 < njmserr) {
				nerr = 13;
				pw->n_err = nerr;
				pw->str_err = "Retrieving diffraction data failed on " + scode + " (" + format("%d", njmserr) + ")";
				goto _exit;
			}
		}
	}

_exit: // final state
	if (pw->b_cancel) { // set cancelled state
		pw->n_state = MS_WS_CANCELLED;
	}
	else {
		if (0 < nerr) { // set finished state
			pw->n_state = MS_WS_ERROR;
		}
		else {
			pw->n_state = MS_WS_FINISHED;
		}
	}
	pw->b_active = false;
	return nerr;
}



unsigned int __cdecl worker_multislice(int gpu_id, int cpu_id, prm_main* pprm, CJMultiSlice* pjms, jms_calc_queue* pq)
{
	unsigned int nerr = 0, ierr = 0; // routine error code
	int whichcode = 0;
	jmsworker w; // current worker data
	size_t itask = 0; // task id
	size_t num_ann = 0; // number of integrating annular detectors
	size_t num_det_slc = 0; // number of detection planes in z
	size_t num_scan = 0; // total number of scan points
	size_t sz_idx = 0; // buffer index
	size_t sz_mem_ann = 0; // memory sizes allocated locally
	size_t i = 0, j = 0; // indices
	unsigned int nyqx = 0, nyqy = 0; // grid nyquist numbers
	float weight = 0.f; // result weight
	float* det_annular = NULL; // pointer to result array for integrating detectors
	float* det_pix_dif = NULL; // pointer to result array for diffraction patterns
	float* det_pix_img = NULL; // pointer to result array for image patterns
	float* dst = NULL;
	float* src = NULL;
	std::string str_pix = "";
	if (NULL == pprm) return 1;
	if (NULL == pjms) return 2;
	if (NULL == pq) return 3;
	num_ann = pprm->detector.v_annular.size();
	num_det_slc = (size_t)pjms->GetDetetionSliceNum();
	num_scan = (size_t)pprm->scan.nx * pprm->scan.ny;
	nyqx = pprm->sample.grid_nx >> 1;
	nyqy = pprm->sample.grid_ny >> 1;
	if (cpu_id >= 0) whichcode = _JMS_CODE_CPU;
	if (gpu_id >= 0) whichcode = _JMS_CODE_GPU;
	if (0 == whichcode) return 4;
	// prepare local result arrays
	if (pprm->detector.b_annular && num_ann > 0 && num_det_slc > 0) {
		// local buffer: holds data of all detectors and all detection planes for one scan position
		sz_mem_ann = sizeof(float) * num_det_slc * num_ann;
		det_annular = (float*)malloc(sz_mem_ann);
		if (NULL == det_annular) {
			nerr = 10;
			goto _exit;
		}
		memset(det_annular, 0, sz_mem_ann);
	}
	if (pprm->detector.b_difpat) {
		det_pix_dif = (float*)pprm->stem_pix_dif.pdata;
	}
	if (pprm->detector.b_image) {
		det_pix_img = (float*)pprm->stem_pix_img.pdata;
	}
	pjms->ResetAveraging(whichcode, cpu_id);
	while (!pq->empty()) {
		itask = pq->request();
		ierr = pq->get_task(itask, w);
		if (0 == ierr) {
			w.whichcode = whichcode;
			w.n_gpu = gpu_id;
			w.n_thread = cpu_id;
			w.pf_res_int = det_annular;
			w.pf_res_dif = det_pix_dif;
			w.pf_res_img = det_pix_img;
			w.id_thread = std::this_thread::get_id();
			pq->set_task_calculate(itask);
			ierr = run_multislice(&w); // run the multislice
			if (ierr == 0 && w.n_state!=MS_WS_CANCELLED) {
				//
				// results per scan point
				// - STEM images of integrated detectors
				if (0 < (w.n_det_flags & (unsigned int)_JMS_DETECT_INTEGRATED)) {
					// write data to result buffer
					// as image series over thickness for each detector
					// - {{plane 1, detector 1},{plane 2, detector 1}, ...{plane N, detector M}}
					for (j = 0; j < num_ann; j++) { // for each detector j
						for (i = 0; i < num_det_slc; i++) { // for each plane i
							sz_idx = (size_t)w.n_scan_idx + (size_t)num_scan * i + (size_t)num_det_slc * num_scan * j;
							dst = &((float*)pprm->stem_images.pdata)[sz_idx];
							src = &w.pf_res_int[i * num_ann + j];
							*dst = *src / w.f_wgt; // write result to pixel
						}
					}
				}
				// - STEM images of integrated detectors
				if (0 < (w.n_det_flags & (unsigned int)_JMS_DETECT_DIFFRACTION)) {
					// write data to result buffers
					// as image series over thickness
					// - position dependent data -> push to disk
					if (pprm->detector.b_difpat) {
						str_pix = format("_px%03d_py%03d", w.n_scan_ix, w.n_scan_iy); // format pixel string
						pprm->stem_pix_dif.shift_org(nyqx + 1, nyqy + 1, 1.f / w.f_wgt); // shift and normalize pattern
						pprm->stem_pix_dif.save(0, "_dif" + str_pix, "bin"); // save pattern
					}
				}
				if (0 < (w.n_det_flags & (unsigned int)_JMS_DETECT_IMAGE)) {
					// write data to result buffers
					// as image series over thickness
					// - position dependent data -> push to disk
					if (pprm->detector.b_image) {
						str_pix = format("_px%03d_py%03d", w.n_scan_ix, w.n_scan_iy); // format pixel string
						pprm->stem_pix_img.normalize(w.f_wgt / pprm->stem_pix_img.f_calc_scale); // normalize image
						pprm->stem_pix_img.save(0, "_img" + str_pix, "bin"); // save image
					}
				}
				//
				pq->set_task_solved(itask); // solved
			}
			else {
				pq->set_task_error(itask); // error
			}
		}
	}
	//
	// collect averaged data from jms
	if (pprm->detector.b_image_avg) { // collect averaged probe image
		size_t sz = pprm->stem_pix_paimg.get_data_bytes();
		float* buf = (float*)malloc(sz);
		ierr = pjms->GetAvgResult(whichcode, _JMS_DETECT_IMAGE_AVG, buf, weight, cpu_id);
		if (0 < ierr) {
			nerr = 60;
			goto _exit;
		}
		ierr = pprm->stem_pix_paimg.add_buffer(buf, weight);
		if (0 < ierr) {
			nerr = 61;
			goto _exit;
		}
		free(buf);
	}
	if (pprm->detector.b_difpat_avg) { // collect averaged diffraction pattern
		size_t sz = pprm->stem_pix_padif.get_data_bytes();
		float* buf = (float*)malloc(sz);
		ierr = pjms->GetAvgResult(whichcode, _JMS_DETECT_DIFFR_AVG, buf, weight, cpu_id);
		if (0 < ierr) {
			nerr = 80;
			goto _exit;
		}
		ierr = pprm->stem_pix_padif.add_buffer(buf, weight);
		if (0 < ierr) {
			nerr = 81;
			goto _exit;
		}
		free(buf);
	}
	//
_exit:
	if (0 < sz_mem_ann) { free(det_annular); det_annular = NULL; sz_mem_ann = 0; }
	return nerr;
}



unsigned int __cdecl singlethread_stem(prm_main *pprm)
{
	if (NULL == pprm) {
		std::cerr << "Error (singlethread_stem): invalid parameter interface." << std::endl;
		return 666;
	}
	unsigned int nerr = 0, ierr = 0;
	int njmserr = 0;
	CJMultiSlice jms;
	jmsworker w;
	unsigned int num_scan = 1;
	unsigned int ix = 0, iy = 0, idx = 0; // scan indices
	unsigned int ix0 = 0, ix1 = pprm->scan.nx, iy0 = 0, iy1 = pprm->scan.ny; // scan range
	unsigned int i = 0, j = 0;
	unsigned int nyqx = pprm->sample.grid_nx >> 1, nyqy = pprm->sample.grid_ny >> 1; // image nyquist numbers
	float scan_rot_sin = 0.f, scan_rot_cos = 1.f;
	float scan_step_x = pprm->scan.get_step_x();
	float scan_step_y = pprm->scan.get_step_y();
	float *det_annular = NULL; // pointer to result array for integrating detectors
	float *det_pix_dif = NULL; // pointer to result array for diffraction patterns
	float *det_pix_img = NULL; // pointer to result array for image patterns
	float *dst = NULL;
	float *src = NULL;
	float weight = 0.f;
	size_t sz_mem_ann = 0; // memory sizes allocated locally
	size_t num_ann = pprm->detector.v_annular.size();
	size_t num_det_slc = 0; // number of detection planes in z
	size_t sz_idx = 0;
	std::string str_pix;
	long long cl_0 = 0, cl_1 = 0;

	if (pprm->btalk) {
		std::cout << std::endl;
		std::cout << "  Initializing single-thread STEM simulation ..." << std::endl;
	}

	cl_0 = pprm->clock.getmsec();

	w.pjms = &jms;

	jms.SetHighTension(pprm->probe.ekv);
	jms.SetGridSize(pprm->sample.grid_nx, pprm->sample.grid_ny);
	jms.SetSupercellSize(pprm->sample.grid_a0, pprm->sample.grid_a1, pprm->sample.grid_a2);
	jms.SetGPUPgrLoading();
	jms.SetPlasmonMC((bool)(pprm->sample.lis_exc_max > 0), pprm->sample.lis_qe, pprm->sample.lis_qc, pprm->sample.lis_mfp, pprm->sample.lis_exc_max);

	if (pprm->gpu_id >= 0) { // run on GPU only
		w.whichcode = _JMS_CODE_GPU;
		njmserr = jms.SetCurrentGPU(pprm->gpu_id);
		if (0 < njmserr) {
			nerr = 1;
			std::cerr << "Error: (singlethread_stem) failed to set GPU #" << pprm->gpu_id << ". (" << njmserr << ")" << std::endl;
			goto _exit;
		}
		w.n_thread = -1;
	}
	else { // run on one CPU process only
		pprm->cpu_num = 1;
		w.whichcode = _JMS_CODE_CPU;
		w.n_thread = 0;
	}

	if (pprm->scan.nx > 0 && pprm->scan.ny > 0) {
		num_scan = pprm->scan.nx * pprm->scan.ny;
		scan_rot_sin = sin(pprm->scan.rotation);
		scan_rot_cos = cos(pprm->scan.rotation);
	}

	// prepare transmission functions
	ierr = prepare_slices(pprm, &jms);
	if (0 < ierr) {
		nerr = 1000 + ierr;
		std::cerr << "Error: (singlethread_stem) failed to prepare transmission functions (" << ierr << ")." << std::endl;
		goto _exit;
	}

	// prepare detectors
	ierr = prepare_detectors(pprm, &jms);
	if (0 < ierr) {
		nerr = 2000 + ierr;
		std::cerr << "Error: (singlethread_stem) failed to prepare detectors (" << ierr << ")." << std::endl;
		goto _exit;
	}
	num_det_slc = (unsigned int)jms.GetDetetionSliceNum();
	w.n_det_flags = pprm->detector.get_jms_flags();

	// prepare cores
	ierr = prepare_cores(pprm, &jms);
	if (0 < ierr) {
		nerr = 3000 + ierr;
		std::cerr << "Error: (singlethread_stem) failed to calculation cores (" << ierr << ")." << std::endl;
		goto _exit;
	}
	
	// prepare probe
	ierr = prepare_probe(pprm, &jms);
	if (0 < ierr) {
		nerr = 4000 + ierr;
		std::cerr << "Error: (singlethread_stem) failed to prepare probe (" << ierr << ")." << std::endl;
		goto _exit;
	}

	// prepare result containers of the prm_main object linked by the function interface
	ierr = pprm->prepare_result_params();
	if (0 < ierr) {
		nerr = 5000 + ierr;
		std::cerr << "Error: (singlethread_stem) failed to prepare result containers." << std::endl;
		goto _exit;
	}

	// prepare local result arrays
	if (pprm->detector.b_annular && num_ann > 0 && num_det_slc > 0) {
		// local buffer: holds data of all detectors and all detection planes for one scan position
		sz_mem_ann = sizeof(float) * num_det_slc * num_ann;
		det_annular = (float*)malloc(sz_mem_ann);
		if (NULL == det_annular) {
			nerr = 10;
			std::cerr << "Error: (singlethread_stem) failed to allocate local scan result buffer." << std::endl;
			goto _exit;
		}
		memset(det_annular, 0, sz_mem_ann);
	}
	if (pprm->detector.b_difpat) {
		det_pix_dif = (float*)pprm->stem_pix_dif.pdata;
	}
	if (pprm->detector.b_image) {
		det_pix_img = (float*)pprm->stem_pix_img.pdata;
	}

	w.n_acc_data = (unsigned int)_JMS_ACCMODE_NONE;
	w.n_repeat = pprm->scan.num_repeat;
	w.b_foc_avg = false;
	w.b_wave_offset = true;
	jms.ResetAveraging(w.whichcode, 0);

	// run the calculation by looping over scan points
	if (pprm->btalk) {
		std::cout << std::endl;
	}
	//
	for (iy = pprm->scan.sub_iy0; iy <= pprm->scan.sub_iy1; iy++) { // slow loop over scan rows
		//
		for (ix = pprm->scan.sub_ix0; ix <= pprm->scan.sub_ix1; ix++) { // fast loop through each scan row
			//
			idx = ix + iy * pprm->scan.nx; // index of the pixel in the scan image
			//
			if (pprm->btalk) {
				std::cout << "  Scanning  x: " << ix+1 << " / " << pprm->scan.nx << "   y: " << iy+1 << " / " << pprm->scan.ny << "     \r";
			}
			//
			// set physical scan pixel position (grid coordinates)
			w.f_dx = pprm->scan.offset_x + scan_rot_cos * scan_step_x * (float)ix - scan_rot_sin * scan_step_y * (float)iy;
			w.f_dy = pprm->scan.offset_y + scan_rot_cos * scan_step_y * (float)iy + scan_rot_sin * scan_step_x * (float)ix;
			w.f_dz = 0.f;
			//
			// set worker parameters
			w.n_scan_idx = idx;
			w.pf_res_int = det_annular;
			w.pf_res_dif = det_pix_dif;
			w.pf_res_img = det_pix_img;
			//
			// run the multislice calculation for this scan pixel
			// (this may include several passes)
			ierr = run_multislice(&w);
			if (0 < ierr) {
				nerr = 20;
				std::cerr << "Error (JMS): " << w.str_err << std::endl;
				std::cerr << "Error (singlethread_stem): failed to run multislice (" << ierr << ")." << std::endl;
				goto _exit;
			}
			//
			// transfer data to result objects
			// - STEM images of integrated detectors
			if (0 < (w.n_det_flags &(unsigned int)_JMS_DETECT_INTEGRATED)) {
				// write data to result buffer
				// as image series over thickness for each detector
				// - {{plane 1, detector 1},{plane 2, detector 1}, ...{plane N, detector M}}
				for (j = 0; j < num_ann; j++) { // for each detector j
					for (i = 0; i < num_det_slc; i++) { // for each plane i
						sz_idx = (size_t)w.n_scan_idx + (size_t)num_scan * i + (size_t)num_det_slc * num_scan * j;
						dst = &((float*)pprm->stem_images.pdata)[sz_idx];
						src = &w.pf_res_int[i * num_ann + j];
						*dst = *src / w.f_wgt;
					}
				}
			}
			// - STEM images of integrated detectors
			if (0 < (w.n_det_flags &(unsigned int)_JMS_DETECT_DIFFRACTION)) {
				// write data to result buffers
				// as image series over thickness
				// - position dependent data -> push to disk
				if (pprm->detector.b_difpat) {
					str_pix = format("_px%03d_py%03d", ix, iy); // format pixel string
					pprm->stem_pix_dif.shift_org(nyqx + 1, nyqy + 1, 1.f / w.f_wgt); // shift and normalize pattern
					pprm->stem_pix_dif.save(0, "_dif" + str_pix, "bin"); // save pattern
				}
			}
			if (0 < (w.n_det_flags &(unsigned int)_JMS_DETECT_IMAGE)) {
				// write data to result buffers
				// as image series over thickness
				// - position dependent data -> push to disk
				if (pprm->detector.b_image) {
					str_pix = format("_px%03d_py%03d", ix, iy); // format pixel string
					pprm->stem_pix_img.normalize(w.f_wgt / pprm->stem_pix_img.f_calc_scale); // normalize image
					pprm->stem_pix_img.save(0, "_img" + str_pix, "bin"); // save image
				}
			}
		} // scan x
	} // scan y
	//
	// collect averaged data from jms
	if (pprm->detector.b_image_avg) { // collect averaged probe image
		ierr = jms.GetAvgResult(w.whichcode, _JMS_DETECT_IMAGE_AVG, (float*)pprm->stem_pix_paimg.pdata, weight, 0);
		if (0 < ierr) {
			nerr = 7001;
			std::cerr << "Error (singlethread_stem): failed to retrieve averaged probe image (" << ierr << ")." << std::endl;
			goto _exit;
		}
		pprm->stem_pix_paimg.f_calc_weight += weight;
	}
	if (pprm->detector.b_difpat_avg) { // collect averaged diffraction pattern
		ierr = jms.GetAvgResult(w.whichcode, _JMS_DETECT_DIFFR_AVG, (float*)pprm->stem_pix_padif.pdata, weight, 0);
		if (0 < ierr) {
			nerr = 7012;
			std::cerr << "Error (singlethread_stem): failed to retrieve averaged probe diffraction (" << ierr << ")." << std::endl;
			goto _exit;
		}
		pprm->stem_pix_padif.f_calc_weight += weight;
	}
	//
	cl_1 = pprm->clock.getmsec();
	if (pprm->btalk) {
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "  Multislice calculation finished in " << std::setprecision(6) << 0.001f * (cl_1 - cl_0) << " s." << std::endl;
	}

_exit:
	if (0 < sz_mem_ann) { free(det_annular); det_annular = NULL; sz_mem_ann = 0; }
	return nerr;
}


unsigned int __cdecl multithread_stem(prm_main *pprm)
{
	if (NULL == pprm) {
		std::cerr << "Error (multithread_stem): invalid parameter interface." << std::endl;
		return 666;
	}

	unsigned int nerr = 0, ierr = 0;
	int njmserr = 0, whichcode = 0, uoff = -1, ncpu = 0;
	unsigned int num_scan = 1;
	unsigned int num_finished = 0, num_failed = 0;
	unsigned int ix = 0, iy = 0, idx = 0; // scan indices
	unsigned int ix0 = 0, ix1 = pprm->scan.nx, iy0 = 0, iy1 = pprm->scan.ny; // scan range
	unsigned int i = 0, j = 0;
	unsigned int nyqx = pprm->sample.grid_nx >> 1, nyqy = pprm->sample.grid_ny >> 1; // image nyquist numbers
	float scan_rot_sin = 0.f, scan_rot_cos = 1.f;
	float scan_step_x = pprm->scan.get_step_x();
	float scan_step_y = pprm->scan.get_step_y();
	float weight = 0.f;
	float perc_solved = 0.f, perc_solved_delta = 0.1f, perc_solved_prev = 0.f;
	size_t num_ann = pprm->detector.v_annular.size();
	size_t num_det_slc = 0; // number of detection planes in z
	size_t sz_idx = 0;
	std::string str_pix;
	long long cl_0 = 0, cl_1 = 0;
	CJMultiSlice jms;
	jms_calc_queue q;
	jmsworker w; // jms worker stats template
	std::thread* pthread = NULL;

	cl_0 = pprm->clock.getmsec();
	
	if (pprm->btalk) {
		std::cout << std::endl;
		std::cout << "  Initializing multi-thread STEM simulation ..." << std::endl;
	}

	// prepare task template
	w.pjms = &jms;
	jms.SetHighTension(pprm->probe.ekv);
	jms.SetGridSize(pprm->sample.grid_nx, pprm->sample.grid_ny);
	jms.SetSupercellSize(pprm->sample.grid_a0, pprm->sample.grid_a1, pprm->sample.grid_a2);
	jms.SetGPUPgrLoading();
	jms.SetPlasmonMC((bool)(pprm->sample.lis_exc_max > 0), pprm->sample.lis_qe, pprm->sample.lis_qc, pprm->sample.lis_mfp, pprm->sample.lis_exc_max);

	whichcode = 0;
	w.n_thread = -1;
	if (pprm->gpu_id >= 0) { // run on GPU only
		whichcode += _JMS_CODE_GPU;
		njmserr = jms.SetCurrentGPU(pprm->gpu_id);
		if (0 < njmserr) {
			nerr = 1;
			std::cerr << "Error (multithread_stem): failed to set GPU #" << pprm->gpu_id << ". (" << njmserr << ")" << std::endl;
			goto _exit;
		}
	}
	ncpu = std::max(0, pprm->cpu_num);
	if (pprm->gpu_id < 0 && ncpu == 0) { 
		ncpu = 1; // fall back to one cpu thread
	}
	if (ncpu > 0) {
		whichcode += _JMS_CODE_CPU;
	}

	if (pprm->scan.nx > 0 && pprm->scan.ny > 0) {
		num_scan = (pprm->scan.sub_ix1 - pprm->scan.sub_ix0 + 1) * (pprm->scan.sub_iy1 - pprm->scan.sub_iy0 + 1);
		scan_rot_sin = sin(pprm->scan.rotation);
		scan_rot_cos = cos(pprm->scan.rotation);
	}
	else {
		nerr = 2;
		std::cerr << "Error (multithread_stem): invalid number of active scan points (" << (pprm->scan.sub_ix1 - pprm->scan.sub_ix0 + 1) << ", " << (pprm->scan.sub_iy1 - pprm->scan.sub_iy0 + 1) << ")." << std::endl;
		goto _exit;
	}

	// prepare transmission functions
	ierr = prepare_slices(pprm, &jms);
	if (0 < ierr) {
		nerr = 1000 + ierr;
		std::cerr << "Error (multithread_stem): failed to prepare transmission functions (" << ierr << ")." << std::endl;
		goto _exit;
	}

	// prepare detectors
	ierr = prepare_detectors(pprm, &jms);
	if (0 < ierr) {
		nerr = 2000 + ierr;
		std::cerr << "Error (multithread_stem): failed to prepare detectors (" << ierr << ")." << std::endl;
		goto _exit;
	}
	num_det_slc = (unsigned int)jms.GetDetetionSliceNum();
	w.n_det_flags = pprm->detector.get_jms_flags();

	// prepare cores
	ierr = prepare_cores(pprm, &jms);
	if (0 < ierr) {
		nerr = 3000 + ierr;
		std::cerr << "Error (multithread_stem): failed to calculation cores (" << ierr << ")." << std::endl;
		goto _exit;
	}

	// prepare probe
	ierr = prepare_probe(pprm, &jms);
	if (0 < ierr) {
		nerr = 4000 + ierr;
		std::cerr << "Error (multithread_stem): failed to prepare probe (" << ierr << ")." << std::endl;
		goto _exit;
	}

	// prepare result containers of the prm_main object linked by the function interface
	ierr = pprm->prepare_result_params();
	if (0 < ierr) {
		nerr = 5000 + ierr;
		std::cerr << "Error (multithread_stem): failed to prepare result containers." << std::endl;
		goto _exit;
	}

	w.n_acc_data = (unsigned int)_JMS_ACCMODE_NONE;
	w.n_repeat = pprm->scan.num_repeat;
	w.b_foc_avg = false;
	w.b_wave_offset = true;

	// initialize queue
	ierr = q.init((size_t)num_scan);
	if (0 < ierr) {
		nerr = 6000 + ierr;
		std::cerr << "Error (multithread_stem): failed to initialize task queue." << std::endl;
		goto _exit;
	}
	// prepare the tasks
	for (iy = pprm->scan.sub_iy0; iy <= pprm->scan.sub_iy1; iy++) { // slow loop over scan rows
		//
		w.n_scan_iy = iy;
		//
		for (ix = pprm->scan.sub_ix0; ix <= pprm->scan.sub_ix1; ix++) { // fast loop through each scan row
			//
			idx = ix + iy * pprm->scan.nx; // index of the pixel in the scan image
			//
			// set physical scan pixel position (grid coordinates)
			w.f_dx = pprm->scan.offset_x + scan_rot_cos * scan_step_x * (float)ix - scan_rot_sin * scan_step_y * (float)iy;
			w.f_dy = pprm->scan.offset_y + scan_rot_cos * scan_step_y * (float)iy + scan_rot_sin * scan_step_x * (float)ix;
			w.f_dz = 0.f;
			//
			// set indices
			w.n_scan_idx = idx;
			w.n_scan_ix = ix;
			//
			// copy to queue
			ierr = q.set_task(idx, w);
			if (0 < ierr) {
				nerr = 6100 + ierr;
				std::cerr << "Error (multithread_stem): failed to set task data in queue at index (" << idx << ")." << std::endl;
				goto _exit;
			}
			//
		} // scan x
	} // scan y
	//
	// start the calculation workers (threads)
	if (pprm->gpu_id >= 0) { // start a gpu thread
		pthread = new std::thread(worker_multislice, pprm->gpu_id, uoff, pprm, &jms, &q);
	}
	if (ncpu > 0) { // start cpu threads
		for (int ithread = 0; ithread < ncpu; ithread++) {
			pthread = new std::thread(worker_multislice, uoff, ithread, pprm, &jms, &q);
		}
	}
	//
	// manage the thread queue progress
	if (pprm->btalk) {
		std::cout << std::endl;
		std::cout << "  Multislice progress: " << format("%5.1f", 0.0) << "%     \r";
	}
	num_finished = (unsigned int)q.solved() + (unsigned int)q.failed();
	while (num_finished < num_scan) { // still calculations running, threads running asynchonously
		if (pprm->btalk) {
			perc_solved = 100.f * (float)num_finished / num_scan;
			if (perc_solved - perc_solved_prev >= perc_solved_delta) {
				std::cout << "  Multislice progress: " << format("%5.1f", perc_solved) << "%     \r";
				perc_solved_prev = perc_solved;
			}
		}
		std::this_thread::sleep_for(std::chrono::milliseconds(5));
		num_finished = (unsigned int)q.solved() + (unsigned int)q.failed();
	}
	if (pprm->btalk) {
		std::cout << "  Multislice progress: " << format("%5.1f", 100.0) << "%       " << std::endl;
	}
	if (q.failed() > 0) {
		nerr = 7061;
		std::cerr << "Error (multithread_stem): " << q.failed() << " of " << q.total() << " calculation tasks failed." << std::endl;
		goto _exit;
	}
	
	//
	cl_1 = pprm->clock.getmsec();
	if (pprm->btalk) {
		std::cout << std::endl;
		std::cout << "  Multislice calculation finished in " << std::setprecision(6) << 0.001f * (cl_1 - cl_0) << " s." << std::endl;
	}

_exit:
	return nerr;
}