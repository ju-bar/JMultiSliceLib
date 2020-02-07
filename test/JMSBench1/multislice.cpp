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


void __cdecl jms_worker_init(jmsworker *pworker)
{
	if (pworker) {
		pworker->b_active = false;
		pworker->b_cancel = false;
		pworker->b_foc_avg = false;
		pworker->b_wave_offset = true;
		pworker->f_dx = 0.f;
		pworker->f_dy = 0.f;
		pworker->f_dz = 0.f;
		pworker->f_fkw = 2.f;
		pworker->f_fs = 0.f;
		pworker->f_wgt = 0.f;
		pworker->id_thread = std::thread::id();
		pworker->n_acc_data = _JMS_ACCMODE_NONE;
		pworker->n_det_flags = _JMS_DETECT_NONE;
		pworker->n_err = 0;
		pworker->n_fkern_num = 7;
		pworker->n_gpu = -1;
		pworker->n_repeat = 1;
		pworker->n_scan_idx = 0;
		pworker->n_state = WS_IDLE;
		pworker->n_thread = -1;
		pworker->pf_res_dif = NULL;
		pworker->pf_res_img = NULL;
		pworker->pf_res_int = NULL;
	}
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
	num_pix = (size_t)(nx * ny);
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



unsigned int __cdecl run_multislice(void* pParam)
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
	jmsworker* pw = NULL;
	if (NULL == pParam) {
		nerr = 1;
		goto _exit;
	}
	// worker state init
	pw = (jmsworker*)pParam;
	pw->b_active = true;
	pw->b_cancel = false;
	pw->n_state = WS_INIT;
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
			pw->n_state = WS_OFFSET_PROBE;
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
		pw->n_state = WS_MULTISLICE;
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
	pw->n_state = WS_READOUT;
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
		pw->n_state = WS_CANCELLED;
	}
	else {
		if (0 < nerr) { // set finished state
			pw->n_state = WS_ERROR;
		}
		else {
			pw->n_state = WS_FINISHED;
		}
	}
	pw->b_active = false;
	return nerr;
}


unsigned int __cdecl singlethread_run(prm_main *pprm)
{
	unsigned int nerr = 0, ierr = 0;
	int njmserr = 0;
	CJMultiSlice jms;
	jmsworker w;
	jms_worker_init(&w);
	unsigned int num_scan = 1;
	unsigned int ix = 0, iy = 0, idx = 0; // scan indices
	unsigned int i = 0, j = 0;
	float scan_rot_sin = 0.f, scan_rot_cos = 1.f;
	float scan_step_x = pprm->scan.get_step_x();
	float scan_step_y = pprm->scan.get_step_y();
	float *det_annular = NULL; // pointer to result array for integrating detectors
	float *dst = NULL;
	float *src = NULL;
	size_t num_mem = 0; // memory size
	size_t num_ann = pprm->detector.v_annular.size();
	size_t num_det_slc = 0; // number of detection planes in z

	if (pprm->btalk) {
		std::cout << std::endl;
		std::cout << "  Initializing single-thread calculation ..." << std::endl;
	}

	w.pjms = &jms;

	jms.SetHighTension(pprm->probe.ekv);
	jms.SetGridSize(pprm->sample.grid_nx, pprm->sample.grid_ny);
	jms.SetGPUPgrLoading();
	// jms.SetPlasmonMC(...);

	if (pprm->gpu_id >= 0) {
		w.whichcode = _JMS_CODE_GPU;
		njmserr = jms.SetCurrentGPU(pprm->gpu_id);
		if (0 < njmserr) {
			nerr = 1;
			std::cerr << "Error: (singlethread_run) failed to set GPU #" << pprm->gpu_id << ". (" << njmserr << ")" << std::endl;
			goto _exit;
		}
	}
	else {
		w.whichcode = _JMS_CODE_CPU;
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
		std::cerr << "Error: (singlethread_run) failed to prepare transmission functions (" << ierr << ")." << std::endl;
		goto _exit;
	}

	// prepare detectors
	ierr = prepare_detectors(pprm, &jms);
	if (0 < ierr) {
		nerr = 2000 + ierr;
		std::cerr << "Error: (singlethread_run) failed to prepare detectors (" << ierr << ")." << std::endl;
		goto _exit;
	}
	num_det_slc = (unsigned int)jms.GetDetetionSliceNum();
	w.n_det_flags = pprm->detector.get_jms_flags();

	// prepare cores
	ierr = prepare_cores(pprm, &jms);
	if (0 < ierr) {
		nerr = 3000 + ierr;
		std::cerr << "Error: (singlethread_run) failed to calculation cores (" << ierr << ")." << std::endl;
		goto _exit;
	}
	
	// prepare probe
	ierr = prepare_probe(pprm, &jms);
	if (0 < ierr) {
		nerr = 4000 + ierr;
		std::cerr << "Error: (singlethread_run) failed to prepare probe (" << ierr << ")." << std::endl;
		goto _exit;
	}

	// prepare local result arrays
	if (pprm->detector.b_annular && num_ann > 0 && num_det_slc > 0) {
		// local buffer: holds data of all detectors and all detection planes for one scan position
		num_mem = sizeof(float) * num_det_slc * num_ann;
		det_annular = (float*)malloc(num_mem);
		if (NULL == det_annular) {
			nerr = 10;
			std::cerr << "Error: (singlethread_run) failed to allocate memory." << std::endl;
			goto _exit;
		}
		memset(det_annular, 0, num_mem);
	}

	// prepare result containers of the prm_main object linked by the function interface
	ierr = pprm->prepare_result_params();
	if (0 < ierr) {
		nerr = 13;
		std::cerr << "Error: (singlethread_run) failed to prepare result containers." << std::endl;
		goto _exit;
	}

	w.n_acc_data = (unsigned int)_JMS_ACCMODE_NONE;
	w.n_repeat = pprm->scan.num_repeat;
	w.b_foc_avg = false;
	w.b_wave_offset = true;

	// run the calculation by looping over scan points
	if (pprm->btalk) {
		std::cout << std::endl;
	}
	//
	for (iy = 0; iy < pprm->scan.ny; iy++) {
		//
		for (ix = 0; ix < pprm->scan.nx; ix++) {
			//
			idx = ix + iy * pprm->scan.nx;
			//
			if (pprm->btalk) {
				std::cout << "  scanning  x: " << ix << " / " << pprm->scan.nx << "   y: " << iy << " / " << pprm->scan.ny << "\r";
			}
			//
			w.f_dx = pprm->scan.offset_x + scan_rot_cos * scan_step_x * (float)ix - scan_rot_sin * scan_step_y * (float)iy;
			w.f_dy = pprm->scan.offset_y + scan_rot_cos * scan_step_y * (float)iy + scan_rot_sin * scan_step_x * (float)ix;
			w.f_dz = 0.f;
			// add random variations here for explicit averaging
			w.n_scan_idx = idx;
			w.pf_res_int = det_annular;
			//
			ierr = run_multislice(&w);
			if (0 < ierr) {
				nerr = 20;
				std::cerr << "Error: (singlethread_run) failed to run multislice (" << ierr << ")." << std::endl;
				goto _exit;
			}
			//
			// transfer results
			if (0 < (w.n_det_flags &(unsigned int)_JMS_DETECT_INTEGRATED)) {
				// write data to result buffer
				// as image series over thickness for each detector
				// - {{plane 1, detector 1},{plane 2, detector 1}, ...{plane N, detector M}}
				for (j = 0; j < num_ann; j++) { // for each detector j
					for (i = 0; i < num_det_slc; i++) { // for each plane i
						dst = &pprm->stem_images.pdata[w.n_scan_idx + i * num_scan + j * num_det_slc * num_scan];
						src = &w.pf_res_int[i + j * num_ann];
						*dst = *src / w.f_wgt;
					}
				}
			}
		}
	}
	if (pprm->btalk) {
		std::cout << std::endl;
	}

_exit:
	if (NULL != det_annular) { free(det_annular); det_annular = NULL; }
	return nerr;
}


unsigned int __cdecl multithread_run(prm_main *pprm)
{
	unsigned int nerr = 0;
	std::cout << "Multi-threaded calculation is not yet implemented." << std::endl;
_exit:
	return nerr;
}