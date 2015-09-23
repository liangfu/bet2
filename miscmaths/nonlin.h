/*    Copyright (C) 2012 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 5.0 (c) 2012, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/9564. */
// Declarations for nonlinear optimisation

#ifndef nonlin_h
#define nonlin_h

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "bfmatrix.h"
#include "newmat.h"

namespace MISCMATHS {

enum NLMethod {NL_VM,                               // Variable-Metric (see NRinC)
               NL_CG,                               // Conjugate-Gradient (see NRinC)
               NL_SCG,                              // Scaled Conjugate-Gradient (See Moller 1993).
               NL_LM,                               // Levenberg-Marquardt (see NRinC)
               NL_GD};                              // Gradient Descent

enum LMType {LM_L, LM_LM};                          // Levenberg or Levenberg-Marquardt

enum VMUpdateType {VM_DFP, VM_BFGS};                // See NRinC chapter 10.

enum CGUpdateType {CG_FR, CG_PR};                   // Fletcher-Reeves, Polak-Ribiere

enum VMMatrixType {VM_OPT,                          // 
                   VM_COL,                          // Store all rank-one updates as column-vectors
                   VM_FULL};                        // Store full estimate of inverse Hessian

enum LinOut {LM_MAXITER,                            // Too many iterations in line-minimisation
             LM_LAMBDA_NILL,                        // Could not find a minima along this direction
             LM_CONV};                              // Line-minimisation converged.

enum NonlinOut {NL_UNDEFINED,                       // Initial value before minimisation
                NL_MAXITER,                         // Too many iterations
                NL_LM_MAXITER,                      // To many iterations during a line-minimisation
                NL_PARCONV,                         // Convergence. Step in parameter space small
                NL_GRADCONV,                        // Convergence. Gradient small
                NL_CFCONV,                          // Convergence. Change in cost-function small
                NL_LCONV};                          // Convergence, lambda very large

const double EPS = 2.0e-16;                         // Losely based on NRinC 20.1

class NonlinException: public std::exception
{
private:
  std::string m_msg;
public:
  NonlinException(const std::string& msg) throw(): m_msg(msg) {}

  virtual const char * what() const throw() {
    return string("Nonlin: msg=" + m_msg).c_str();
  }

  ~NonlinException() throw() {}
};

// NonlinParam is a struct that contains the
// information about "how" the 
// minisation should be performed. I.e. it
// contains things like choice of minimisation
// algorithm, # of parameters, converegence 
// criteria etc

class NonlinParam
{
public:

  NonlinParam(int                   pnpar,
              NLMethod              pmtd,
              NEWMAT::ColumnVector  ppar=NEWMAT::ColumnVector(),
              bool                  plogcf=false,
              bool                  ploglambda=false,
              bool                  plogpar=false, 
              int                   pmaxiter=200,
              double                pcftol=1.0e-8,
              double                pgtol=1.0e-8,
              double                pptol=4.0*EPS,
              VMUpdateType          pvmut=VM_BFGS,
              double                palpha=1.0e-4,
              double                pstepmax=10,
              int                   plm_maxiter=50,
              int                   pmaxrestart=0,
              bool                  pautoscale=true,
              CGUpdateType          pcgut=CG_PR,
              double                plm_ftol=1.0e-3,
              LMType                plmtype=LM_LM,
              double                pltol=1.0e20,
	      int                   pcg_maxiter=200,
	      double                pcg_tol=1.0e-6,
              double                plambda=0.1)
    : npar(pnpar), mtd(pmtd), logcf(plogcf), 
    loglambda(ploglambda), logpar(plogpar), maxiter(pmaxiter), 
    cftol(pcftol), gtol(pgtol), ptol(pptol), vmut(pvmut), 
    alpha(palpha), stepmax(pstepmax), lm_maxiter(plm_maxiter),
    maxrestart(pmaxrestart), autoscale(pautoscale), cgut(pcgut), 
    lm_ftol(plm_ftol), lmtype(plmtype), ltol(pltol), cg_maxiter(pcg_maxiter), 
    cg_tol(pcg_tol), lambda(), cf(), par(), niter(0), nrestart(0), status(NL_UNDEFINED)
  {
    lambda.push_back(plambda);
    if (ppar.Nrows()) SetStartingEstimate(ppar);
    else {
      NEWMAT::ColumnVector  tmp(npar);
      tmp = 0.0;
      SetStartingEstimate(tmp);
    }
  }
  ~NonlinParam() {}                  

  // Routines to check values

  int NPar() const {return(npar);}
  NLMethod Method() const {return(mtd);}
  int MaxIter() const {return(maxiter);}
  int NIter() const {return(niter);}
  double FractionalCFTolerance() const {return(cftol);}
  double FractionalGradientTolerance() const {return(gtol);}
  double FractionalParameterTolerance() const {return(ptol);}
  VMUpdateType VariableMetricUpdate() const {return(vmut);}
  double VariableMetricAlpha() const {return(alpha);}
  int MaxVariableMetricRestarts() const {return(maxrestart);}
  int VariableMetricRestarts() const {return(nrestart);}
  bool VariableMetricAutoScale() const {return(autoscale);}
  double LineSearchMaxStep() const {return(stepmax);}
  int LineSearchMaxIterations() const {return(lm_maxiter);}
  CGUpdateType ConjugateGradientUpdate() const {return(cgut);}
  double LineSearchFractionalParameterTolerance() const {return(lm_ftol);}
  LMType GaussNewtonType() const {return(lmtype);}
  double LambdaConvergenceCriterion() const {return(ltol);}
  int EquationSolverMaxIter() const {return(cg_maxiter);}
  double EquationSolverTol() const {return(cg_tol);}
  bool LoggingParameters() const {return(logpar);}
  bool LoggingCostFunction() const {return(logcf);}
  bool LoggingLambda() const {return(loglambda);}

  // Routines to get output

  double Lambda() const {return(lambda.back());}
  double InitialLambda() const {if (loglambda) return(lambda[0]); else {throw NonlinException("InitialLabda: Lambda not logged"); return(0.0);}}
  const std::vector<double>& LambdaHistory() const {if (loglambda) return(lambda); else {throw NonlinException("InitialLabda: Lambda not logged"); return(lambda);}}
  const NEWMAT::ColumnVector& Par() const {return(par.back());}
  const NEWMAT::ColumnVector& InitialPar() const {if (logpar) return(par[0]); else {throw NonlinException("InitialPar: Parameters not logged"); return(par[0]);}}
  const std::vector<NEWMAT::ColumnVector>& ParHistory() const {if (logpar) return(par); else {throw NonlinException("ParHistory: Parameters not logged"); return(par);}}
  double CF() const {return(cf.back());}
  double InitialCF() const {if (logcf) return(cf[0]); else {throw NonlinException("InitialCF: Cost-function not logged"); return(cf[0]);}}
  const std::vector<double> CFHistory() const {if (logcf) return(cf); else {throw NonlinException("CFHistory: Cost-function not logged"); return(cf);}}
  NonlinOut Status() const {return(status);}
  bool Success() const { switch(status) { case NL_UNDEFINED: case NL_MAXITER: case NL_LM_MAXITER: return(false); break; default: return(true); } };
  std::string TextStatus() const;
      
  // Routines to set values of steering parameters
  void SetMethod(NLMethod pmtd) {mtd = pmtd;}
  void LogCF(bool flag=true) {logcf = flag;}
  void LogPar(bool flag=true) {logpar = flag;}
  void LogLambda(bool flag=true) {loglambda = flag;}
  void SetStartingEstimate(NEWMAT::ColumnVector& sp) {
    if (niter) throw NonlinException("SetStartingEstimates: Object has to be reset before setting new starting parameters");
    SetPar(sp);
  }
  void SetMaxIter(unsigned int pmiter) {maxiter = pmiter;}
  void SetFractionalCFTolerance(double pcftol) {
    if (pcftol>0.5) throw NonlinException("SetFractionalCFTolerance: Nonsensically large tolerance");
    else if (pcftol <= 0.0) NonlinException("SetFractionalCFTolerance: Tolerance must be non-zero and positive");
    cftol = pcftol;
  }
  void SetFractionalGradientTolerance(double pgtol) {
    if (pgtol>0.5) throw NonlinException("SetFractionalGradientTolerance: Nonsensically large tolerance");
    else if (pgtol <= 0.0) NonlinException("SetFractionalGradientTolerance: Tolerance must be non-zero and positive");
    gtol = pgtol;
  }
  void SetFractionalParameterTolerance(double pptol) {
    if (pptol>0.5) throw NonlinException("SetFractionalParameterTolerance: Nonsensically large tolerance");
    else if (pptol <= 0.0) NonlinException("SetFractionalParameterTolerance: Tolerance must be non-zero and positive");
    ptol = pptol;
  }
  void SetVariableMetricUpdate(VMUpdateType pvmut) {vmut = pvmut;}
  void SetVariableMetricAlpha(double palpha) {
    if (palpha>=1.0 || palpha<=0.0) throw NonlinException("SetVariableMetricAlpha: Alpha must be between 0 and 1");
    alpha = palpha;
  }
  void SetMaxVariableMetricRestarts(unsigned int pmaxrestart) {maxrestart = pmaxrestart;}
  void SetVariableMetricAutoScale(bool flag=true) {autoscale = flag;}
  void SetLineSearchMaxStep(double pstepmax) {
    if (pstepmax<=0) throw NonlinException("SetLineSearchMaxStep: maxstep must be non-zero and positive");
    stepmax = pstepmax;
  }
  void SetLineMinimisationMaxIterations(unsigned int plm_maxiter) {lm_maxiter = plm_maxiter;}
  void SetConjugateGradientUpdate(CGUpdateType pcgut) {cgut = pcgut;}
  void SetLineMinimisationFractionalParameterTolerance(double plm_ftol) {
    if (plm_ftol>0.5) throw NonlinException("SetLineMinimisationFractionalParameterTolerance: Nonsensically large tolerance");
    else if (plm_ftol <= 0.0) NonlinException("SetLineMinimisationFractionalParameterTolerance: Tolerance must be non-zero and positive");
    lm_ftol = plm_ftol;
  }
  void SetGaussNewtonType(LMType plmtype) {lmtype = plmtype;}
  void SetLambdaConvergenceCriterion(double pltol) {
    if (pltol<1.0) throw NonlinException("SetLambdaConvergenceCriterion: Nonsensically small tolerance");
    ltol = pltol;
  }
  void SetEquationSolverMaxIter(int pcg_maxiter) {cg_maxiter = pcg_maxiter;}
  void SetEquationSolverTol(double pcg_tol) {cg_tol = pcg_tol;}

  
  // Reset is used to reset a NonlinParam object after it has run to convergence, thereby allowing it
  // to be reused with a different CF object. This is to avoid the cost of creating the object many
  // times when fitting for example multiple voxels.
  void Reset() {}
  // Routines used by the (global) non-linear fitting routines. Note that these can
  // all be called for const objects.
  void SetPar(const NEWMAT::ColumnVector& p) const {
    if (p.Nrows() != npar) throw NonlinException("SetPar: Mismatch between starting vector and # of parameters");
    if (logpar || !par.size()) par.push_back(p); 
    else par[0] = p;
  }
  void SetCF(double pcf) const {
    if (logcf || !cf.size()) cf.push_back(pcf); 
    else cf[0] = pcf;
  }
  void SetLambda(double pl) const {
    if (loglambda || !lambda.size()) lambda.push_back(pl); 
    else lambda[0] = pl;
  }
  bool NextIter(bool success=true) const {if (success && niter++ >= maxiter) return(false); else return(true);}
  bool NextRestart() const {if (nrestart++ >= maxrestart) return(false); else return(true);}
  void SetStatus(NonlinOut  pstatus) const {status = pstatus;}
private:

  //         INPUT PARAMETERS
  //
  //         Paramaters that apply to all algorithms

  const int                  npar;       // # of parameters
  NLMethod                   mtd;        // Minimisation method
  bool                       logcf;      // If true, history of cost-function is logged
  bool                       loglambda;  // If true, history of lambda is logged
  bool                       logpar;     // If true history of parameters is logged
  int                        maxiter;    // Maximum # of iterations allowed
  double                     cftol;      // Tolerance for cost-function gonvergence criterion
  double                     gtol;       // Tolerance for gradient convergence criterion
  double                     ptol;       // Tolerance for parameter convergence criterion

  //         Parameters that apply to Variable-Metric Algorithm

  VMUpdateType               vmut;       // DFP or BFGS
  double                     alpha;      // Criterion for convergence in line minimisation
  double                     stepmax;    // Maximum step length for line minimisation
  int                        lm_maxiter; // Maximum # of iterations for line minimisation
  int                        maxrestart; // Maximum # of restarts that should be done.      
  bool                       autoscale;  // "Automatic" search for optimal scaling

  //         Parameters that apply to CG algorithm

  CGUpdateType               cgut;       // Fletcher-Reeves or Polak-Ribiere
  double                     lm_ftol;    // Convergence criterion for line-search

  //         Parameters that apply to LM algorithm

  LMType                     lmtype;     // Levenberg or Levenberg-Marquardt
  double                     ltol;       // Convergence criterion based on large lambda
  int                        cg_maxiter; // Maximum # of iterations for iterative "inverse" of Hessian
  double                     cg_tol;     // Tolerance for iterative "inverse" of Hessian

  //
  //         OUTPUT PARAMETERS
  //
  mutable std::vector<double>                        lambda;     // (History of) lambda (LM and SCG type minimisation)
  mutable std::vector<double>                        cf;         // (History of) cost-function
  mutable std::vector<NEWMAT::ColumnVector>          par;        // (History of) Parameter estimates
  mutable int                                        niter;      // Number of iterations
  mutable int                                        nrestart;   // Number of restarts
  mutable NonlinOut                                  status;     // Output status  
  
  NonlinParam& operator=(const NonlinParam& rhs); // Hide assignment

};
  
// NonlinCF (Cost Function) is a virtual
// class that defines a minimal interface.
// By subclassing NonlinCF the "user" can
// create a class that allows him/her to
// use NONLIN to minimise his/her function.
class NonlinCF
{
private:
  NonlinCF& operator=(const NonlinCF& rhs); // Hide assignment
public:
  NonlinCF() {}
  virtual ~NonlinCF() {}
  virtual double sf() const {return(1.0);}
  virtual NEWMAT::ReturnMatrix grad(const NEWMAT::ColumnVector& p) const;  
  virtual boost::shared_ptr<BFMatrix> hess(const NEWMAT::ColumnVector& p,
                                           boost::shared_ptr<BFMatrix> iptr=boost::shared_ptr<BFMatrix>()) const;  
  virtual double cf(const NEWMAT::ColumnVector& p) const = 0;              
};

// Varmet matrix is a "helper" class
// that makes it a little easier to
// implement variable-metric minimisation.
class VarmetMatrix
{
private:
  int                                 sz;
  VMMatrixType                        mtp;
  VMUpdateType                        utp;
  NEWMAT::Matrix                      mat;
  std::vector<double>                 sf;
  std::vector<NEWMAT::ColumnVector>   vec;

  VarmetMatrix& operator=(const VarmetMatrix& rhs);  // Hide assignment

public:
  explicit VarmetMatrix(int psz, VMMatrixType pmtp, VMUpdateType putp) 
  : sz(psz), mtp(pmtp), utp(putp)
  {
    if (sz > 0 && mtp == VM_OPT) {
      if (sz < 100) {
        mtp = VM_FULL;
	NEWMAT::IdentityMatrix  tmp(sz);
        mat = tmp;
      }
      else {
        mtp = VM_COL;
      }
    }
  }

  ~VarmetMatrix() {}

  int size() {return(sz);}
  VMUpdateType updatetype() {return(utp);}
  VMMatrixType matrixtype() {return(mtp);}
  void print() const;

  void reset()
  {
    if (sz > 0) {
      if (mtp == VM_FULL) {
	NEWMAT::IdentityMatrix  tmp(sz);
        mat = tmp;
      }
      else if (mtp == VM_COL) {
        sf.clear();
        vec.clear();
      }
    }
  }

  void update(const NEWMAT::ColumnVector& pdiff,  // x_{i+1} - x_i
              const NEWMAT::ColumnVector& gdiff); // \nabla f_{i+1} - \nabla f_i
  friend NEWMAT::ColumnVector operator*(const VarmetMatrix&          m,
                                        const NEWMAT::ColumnVector&  v); 
};

// Declaration of (global) main function for minimisation

NonlinOut nonlin(const NonlinParam& p, const NonlinCF& cfo);

// Declaration of global utility functions

pair<NEWMAT::ColumnVector,NEWMAT::ColumnVector> check_grad(const NEWMAT::ColumnVector&  par,
                                                           const NonlinCF&                     cfo);

pair<boost::shared_ptr<BFMatrix>,boost::shared_ptr<BFMatrix> > check_hess(const NEWMAT::ColumnVector& par,
                                                                          const NonlinCF&     cfo);

} // End namespace MISCMATHS

#endif // end #ifndef nonlin_h
