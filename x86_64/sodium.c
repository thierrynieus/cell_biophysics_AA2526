/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__sodium
#define _nrn_initial _nrn_initial__sodium
#define nrn_cur _nrn_cur__sodium
#define _nrn_current _nrn_current__sodium
#define nrn_jacob _nrn_jacob__sodium
#define nrn_state _nrn_state__sodium
#define _net_receive _net_receive__sodium 
#define rates rates__sodium 
#define states states__sodium 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg(int);
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gnabar _p[0]
#define gnabar_columnindex 0
#define Am _p[1]
#define Am_columnindex 1
#define v0_Am _p[2]
#define v0_Am_columnindex 2
#define k_Am _p[3]
#define k_Am_columnindex 3
#define Bm _p[4]
#define Bm_columnindex 4
#define v0_Bm _p[5]
#define v0_Bm_columnindex 5
#define k_Bm _p[6]
#define k_Bm_columnindex 6
#define Ah _p[7]
#define Ah_columnindex 7
#define v0_Ah _p[8]
#define v0_Ah_columnindex 8
#define k_Ah _p[9]
#define k_Ah_columnindex 9
#define Bh _p[10]
#define Bh_columnindex 10
#define v0_Bh _p[11]
#define v0_Bh_columnindex 11
#define k_Bh _p[12]
#define k_Bh_columnindex 12
#define ina _p[13]
#define ina_columnindex 13
#define minf _p[14]
#define minf_columnindex 14
#define hinf _p[15]
#define hinf_columnindex 15
#define mtau _p[16]
#define mtau_columnindex 16
#define htau _p[17]
#define htau_columnindex 17
#define gna _p[18]
#define gna_columnindex 18
#define m _p[19]
#define m_columnindex 19
#define h _p[20]
#define h_columnindex 20
#define Dm _p[21]
#define Dm_columnindex 21
#define Dh _p[22]
#define Dh_columnindex 22
#define q10 _p[23]
#define q10_columnindex 23
#define v _p[24]
#define v_columnindex 24
#define _g _p[25]
#define _g_columnindex 25
#define _ion_ina	*_ppvar[0]._pval
#define _ion_dinadv	*_ppvar[1]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_rates(void);
 static void _hoc_vtrap(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_sodium", _hoc_setdata,
 "rates_sodium", _hoc_rates,
 "vtrap_sodium", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_sodium
 extern double vtrap( _threadargsprotocomma_ double , double );
 /* declare global and static user variables */
#define ena ena_sodium
 double ena = 55;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "ena_sodium", "mV",
 "gnabar_sodium", "mho/cm2",
 "ina_sodium", "mA/cm2",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "ena_sodium", &ena_sodium,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"sodium",
 "gnabar_sodium",
 "Am_sodium",
 "v0_Am_sodium",
 "k_Am_sodium",
 "Bm_sodium",
 "v0_Bm_sodium",
 "k_Bm_sodium",
 "Ah_sodium",
 "v0_Ah_sodium",
 "k_Ah_sodium",
 "Bh_sodium",
 "v0_Bh_sodium",
 "k_Bh_sodium",
 0,
 "ina_sodium",
 "minf_sodium",
 "hinf_sodium",
 "mtau_sodium",
 "htau_sodium",
 "gna_sodium",
 0,
 "m_sodium",
 "h_sodium",
 0,
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 26, _prop);
 	/*initialize range parameters*/
 	gnabar = 0.12;
 	Am = 0.1;
 	v0_Am = 40;
 	k_Am = 10;
 	Bm = 4;
 	v0_Bm = 65;
 	k_Bm = 18;
 	Ah = 0.07;
 	v0_Ah = 65;
 	k_Ah = 20;
 	Bh = 1;
 	v0_Bh = 35;
 	k_Bh = 10;
 	_prop->param = _p;
 	_prop->param_size = 26;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 	_ppvar[0]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[1]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _sodium_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 26, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 sodium /home/tnieus/Documents/teaching/AA 2025-2026/Cell Biophysics/lezioni/cell_biophysics_AA2526/sodium.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "sodium current";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
   rates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0 ) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0 ) ) ) / htau ) - h) ;
   }
  return 0;
}
 
static int  rates ( _threadargsprotocomma_ double _lv ) {
   double _lalpha_m , _lbeta_m , _lalpha_h , _lbeta_h ;
 _lalpha_m = Am * vtrap ( _threadargscomma_ - ( _lv + v0_Am ) , k_Am ) ;
   _lbeta_m = Bm * exp ( - ( _lv + v0_Bm ) / k_Bm ) ;
   mtau = 1.0 / ( q10 * ( _lalpha_m + _lbeta_m ) ) ;
   minf = _lalpha_m / ( _lalpha_m + _lbeta_m ) ;
   _lalpha_h = Ah * exp ( - ( _lv + v0_Ah ) / k_Ah ) ;
   _lbeta_h = Bh / ( 1.0 + exp ( - ( _lv + v0_Bh ) / k_Bh ) ) ;
   htau = 1.0 / ( q10 * ( _lalpha_h + _lbeta_h ) ) ;
   hinf = _lalpha_h / ( _lalpha_h + _lbeta_h ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap ( _threadargsprotocomma_ double _lx , double _ly ) {
   double _lvtrap;
 if ( fabs ( _lx / _ly ) < 1e-6 ) {
     _lvtrap = _ly * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _lvtrap = _lx / ( exp ( _lx / _ly ) - 1.0 ) ;
     }
   
return _lvtrap;
 }
 
static void _hoc_vtrap(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  m = m0;
 {
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
   q10 = pow( 3.0 , ( ( celsius - 6.3 ) / 10.0 ) ) ;
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gna = gnabar * pow( m , 3.0 ) * h ;
   ina = gna * ( v - ena ) ;
   }
 _current += ina;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = m_columnindex;  _dlist1[0] = Dm_columnindex;
 _slist1[1] = h_columnindex;  _dlist1[1] = Dh_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/tnieus/Documents/teaching/AA 2025-2026/Cell Biophysics/lezioni/cell_biophysics_AA2526/sodium.mod";
static const char* nmodl_file_text = 
  "TITLE sodium current\n"
  " \n"
  "COMMENT\n"
  "Based on https://neuron.yale.edu/neuron/docs/hodgkin-huxley-using-rxd\n"
  "\n"
  "\n"
  "q10 = 3^((celsius - 6.3)/10)\n"
  ":\"m\" sodium activation system\n"
  "alpha = .1 * vtrap(-(v+40),10)\n"
  "beta =  4 * exp(-(v+65)/18)\n"
  "sum = alpha + beta\n"
  "mtau = 1/(q10*sum)\n"
  "minf = alpha/sum\n"
  ":\"h\" sodium inactivation system\n"
  "alpha = .07 * exp(-(v+65)/20)\n"
  "beta = 1 / (exp(-(v+35)/10) + 1)\n"
  "sum = alpha + beta\n"
  "htau = 1/(q10*sum)\n"
  "hinf = alpha/sum\n"
  "\n"
  "\n"
  ":\"n\" potassium activation system\n"
  "alpha = .01*vtrap(-(v+55),10) \n"
  "beta = .125*exp(-(v+65)/80)\n"
  "sum = alpha + beta\n"
  "ntau = 1/(q10*sum)\n"
  "ninf = alpha/sum\n"
  "\n"
  "\n"
  "\n"
  "ENDCOMMENT\n"
  " \n"
  "UNITS {\n"
  "    (mA) = (milliamp)\n"
  "    (mV) = (millivolt)	\n"
  "}\n"
  " \n"
  "NEURON {\n"
  " 	SUFFIX sodium\n"
  "	USEION na WRITE ina\n"
  "	RANGE gnabar, gna, minf, hinf, mtau, htau, m, h, ina, Am, v0_Am, k_Am, Bm, v0_Bm, k_Bm, Ah, v0_Ah, k_Ah, Bh, v0_Bh, k_Bh\n"
  "}\n"
  " \n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  " \n"
  "PARAMETER {\n"
  "  	\n"
  "	ena	= 55	(mV)    : CHECK\n"
  "	gnabar	= 0.12 (mho/cm2) : OK\n"
  "	\n"
  "	: HH params\n"
  "	Am = 0.1\n"
  "    v0_Am = 40\n"
  "    k_Am = 10\n"
  "    Bm = 4\n"
  "    v0_Bm = 65\n"
  "    k_Bm = 18 \n"
  "    \n"
  "	Ah = 0.07\n"
  "    v0_Ah = 65\n"
  "    k_Ah = 20\n"
  "    Bh = 1\n"
  "    v0_Bh = 35\n"
  "    k_Bh = 10\n"
  "     \n"
  "}\n"
  " \n"
  "STATE {\n"
  "    m \n"
  "    h\n"
  "}\n"
  " \n"
  "ASSIGNED {\n"
  "    v           (mV)\n"
  "    ina         (mA/cm2)\n"
  "    celsius		(degC)\n"
  " 	minf\n"
  "    hinf\n"
  "    mtau\n"
  "    htau\n"
  "    gna\n"
  "    q10\n"
  "}\n"
  " \n"
  "BREAKPOINT {\n"
  "    SOLVE states METHOD cnexp\n"
  "    gna = gnabar * m^3 * h\n"
  "    ina = gna * (v - ena)\n"
  "}\n"
  " \n"
  "UNITSOFF\n"
  " \n"
  "INITIAL {\n"
  "    rates(v)\n"
  "    m = minf\n"
  "    h = hinf\n"
  "    q10 = 3.0^((celsius - 6.3)/10.0)\n"
  "}\n"
  "\n"
  "DERIVATIVE states {  \n"
  "    :Computes states variable m and h at the current v and t.\n"
  "    rates(v)\n"
  "    m' = (minf - m) / mtau\n"
  "    h' = (hinf - h) / htau\n"
  "}\n"
  " \n"
  "PROCEDURE rates(v) {  \n"
  "    :Computes rate and other constants at current v.\n"
  "    :Call once from HOC to initialize inf at resting v.\n"
  "    LOCAL alpha_m, beta_m, alpha_h, beta_h\n"
  "    : m\n"
  "    alpha_m = Am * vtrap(-(v + v0_Am), k_Am)\n"
  "    beta_m = Bm * exp(-(v + v0_Bm) / k_Bm)\n"
  "    mtau = 1.0/(q10 * (alpha_m + beta_m))\n"
  "    minf = alpha_m / (alpha_m + beta_m)\n"
  "    : h\n"
  "    alpha_h = Ah * exp(-(v + v0_Ah)/k_Ah)\n"
  "    beta_h = Bh / ( 1 + exp(-(v + v0_Bh)/k_Bh))\n"
  "    htau = 1 / ( q10 * ( alpha_h + beta_h ) )\n"
  "    hinf = alpha_h / ( alpha_h + beta_h )    \n"
  "\n"
  "  \n"
  "}\n"
  "\n"
  "FUNCTION vtrap(x (mV),y (mV)) (mV) {\n"
  "    : vtrap(x,y) is 1/(exp(x/y)-1) if |x/y|>=1e-6 or y*(1.0 - x/y/2.0) otherwise.\n"
  "    if (fabs(x/y) < 1e-6) {\n"
  "            vtrap = y*(1 - x/y/2)\n"
  "    }else{\n"
  "            vtrap = x/(exp(x/y) - 1)\n"
  "    }\n"
  "}\n"
  "\n"
  "\n"
  " \n"
  "UNITSON\n"
  ;
#endif
