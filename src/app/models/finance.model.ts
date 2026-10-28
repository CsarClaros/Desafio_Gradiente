export type ScenarioKey = 'ideal' | 'stress' | 'ill';
export type Preconditioner = 'jacobi' | 'ssor' | 'none';

export interface FinancialScenario {
  key: ScenarioKey;
  title: string;
  subtitle: string;
  justification: string;
  financialReading: string;
  matrix: string[][];
  rhs: string[];
}

export interface Iteration {
  k: number;
  x: number[];
  residual: number;
  alpha: number | null;
  beta: number | null;
}

export interface SolverResult {
  x: number[];
  iterations: number;
  converged: boolean;
  residual: number;
  history: Iteration[];
  message: string;
}

export interface ErrorRow {
  label: string;
  reference: number;
  pcg: number;
  abs: number;
  relative: number | null;
  percent: number | null;
}

export interface Analysis {
  valid: boolean;
  message: string;
  symmetric: boolean;
  spd: boolean;
  n: number;
  A: number[][];
  b: number[];
  minors: number[];
  eigenvalues: number[];
  eigenPre: number[];
  cond: number | null;
  condPre: number | null;
  condSym: number | null;
  reference: number[];
  preMatrix: number[][];
  initialResidual: number[];
  cg: SolverResult | null;
  pcg: SolverResult | null;
  errors: ErrorRow[];
  referenceResidual: number | null;
  preconditioner: Preconditioner;
}

export interface QAResult {
  code: string;
  name: string;
  input: string;
  expected: string;
  obtained: string;
  passed: boolean;
  conclusion: string;
}
