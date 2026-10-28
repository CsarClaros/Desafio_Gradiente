import { Injectable } from '@angular/core';
import { Analysis, ErrorRow, FinancialScenario, Iteration, Preconditioner, QAResult, ScenarioKey, SolverResult } from '../models/finance.model';

export type Matrix = number[][];
const EPS = 1e-14;

function zeros(n: number): number[] { return Array(n).fill(0); }
function dot(a: number[], b: number[]): number { return a.reduce((s, v, i) => s + v * b[i], 0); }
function norm(a: number[]): number { return Math.sqrt(dot(a, a)); }
function matVec(a: Matrix, x: number[]): number[] { return a.map(row => dot(row, x)); }
function sub(a: number[], b: number[]): number[] { return a.map((v, i) => v - b[i]); }
function transpose(a: Matrix): Matrix { return a[0].map((_, i) => a.map(row => row[i])); }
function matMul(a: Matrix, b: Matrix): Matrix { const bt = transpose(b); return a.map(row => bt.map(col => dot(row, col))); }
function identity(n: number): Matrix { return Array.from({length:n}, (_,i) => Array.from({length:n}, (_,j) => i===j ? 1 : 0)); }
function isFiniteMatrix(a: Matrix): boolean { return a.every(row => row.every(Number.isFinite)); }

/** Factorización de Cholesky; los menores principales vienen del producto de pivotes L_ii². */
export function cholesky(a: Matrix): {L: Matrix; minors:number[]} | null {
  const n=a.length, L:Matrix=Array.from({length:n},()=>zeros(n));
  const minors:number[]=[]; let det=1;
  for(let i=0;i<n;i++){
    for(let j=0;j<=i;j++){
      let value=a[i][j];
      for(let k=0;k<j;k++) value-=L[i][k]*L[j][k];
      if(i===j){
        if(!Number.isFinite(value) || value<=Math.max(1e-13,Math.abs(a[i][i])*1e-13)) return null;
        L[i][i]=Math.sqrt(value); det*=value; minors.push(det);
      } else { L[i][j]=value/L[j][j]; }
    }
  }
  return {L,minors};
}
function forward(L:Matrix,b:number[]):number[] {
  const x=zeros(b.length);
  for(let i=0;i<b.length;i++){
    let s=b[i]; for(let j=0;j<i;j++)s-=L[i][j]*x[j]; x[i]=s/L[i][i];
  }
  return x;
}
function backUpper(U:Matrix,b:number[]):number[]{
  const n=b.length,x=zeros(n);
  for(let i=n-1;i>=0;i--){let s=b[i]; for(let j=i+1;j<n;j++)s-=U[i][j]*x[j]; x[i]=s/U[i][i];}
  return x;
}
function solveCholeskyFactor(L:Matrix,b:number[]):number[]{ return backUpper(transpose(L),forward(L,b)); }
function inverseLower(L:Matrix):Matrix {
  const n=L.length, inv=Array.from({length:n},()=>zeros(n));
  for(let j=0;j<n;j++){
    const e=zeros(n);e[j]=1;
    const col=forward(L,e);
    for(let i=0;i<n;i++) inv[i][j]=col[i];
  }
  return inv;
}
/** Autovalores de una matriz simétrica vía rotaciones ortogonales de Jacobi. */
export function eigenSym(a:Matrix):number[]{
  const n=a.length, m=a.map(row=>row.slice());
  const scale=Math.max(1,...a.map(row=>Math.max(...row.map(Math.abs))));
  for(let it=0;it<100*n*n;it++){
    let p=0,q=1,max=0;
    for(let i=0;i<n;i++) for(let j=i+1;j<n;j++){
      const v=Math.abs(m[i][j]);if(v>max){max=v;p=i;q=j;}
    }
    if(max<=scale*2e-15) break;
    const tau=(m[q][q]-m[p][p])/(2*m[p][q]);
    const t=(tau>=0?1:-1)/(Math.abs(tau)+Math.sqrt(1+tau*tau));
    const c=1/Math.sqrt(1+t*t),s=t*c;
    const app=m[p][p],aqq=m[q][q],apq=m[p][q];
    m[p][p]=app-t*apq; m[q][q]=aqq+t*apq;
    m[p][q]=0;m[q][p]=0;
    for(let k=0;k<n;k++) if(k!==p && k!==q){
      const akp=m[k][p],akq=m[k][q];
      m[k][p]=c*akp-s*akq; m[p][k]=m[k][p];
      m[k][q]=s*akp+c*akq; m[q][k]=m[k][q];
    }
  }
  return m.map((row,i)=>row[i]).sort((x,y)=>x-y);
}
function precondMatrix(A:Matrix, method:Preconditioner,omega:number):Matrix {
  const n=A.length;
  if(method==='none')return identity(n);
  if(method==='jacobi')return A.map((row,i)=>row.map((_,j)=>i===j?A[i][i]:0));
  // M_SSOR = (D+ωL) D⁻¹ (D+ωU) / [ω(2−ω)], 0<ω<2.
  const left=A.map((row,i)=>row.map((v,j)=>j<i?omega*v:j===i?v:0));
  const right=transpose(left);
  const mid=left.map((row,i)=>row.map((v,j)=>v/A[j][j]));
  return matMul(mid,right).map(row=>row.map(v=>v/(omega*(2-omega))));
}
function conditionSym(a:Matrix):{cond:number,eig:number[]}{
  const eig=eigenSym(a), min=eig[0],max=eig[eig.length-1];
  return {cond:min>0?max/min:Number.POSITIVE_INFINITY,eig};
}
function relativeResidual(A:Matrix,b:number[],x:number[]):number {
  const denom=norm(b); const actual=norm(sub(b,matVec(A,x)));
  return denom>0?actual/denom:actual;
}
function runCG(A:Matrix,b:number[],x0:number[],tol:number,maxIter:number,preL:Matrix|null):SolverResult{
  const n=A.length; let x=x0.slice(),r=sub(b,matVec(A,x));
  const apply=(v:number[]):number[]=>preL?solveCholeskyFactor(preL,v):v.slice();
  let z=apply(r),p=z.slice(),rz=dot(r,z);
  let rr=relativeResidual(A,b,x);
  const history:Iteration[]=[{k:0,x:x.slice(),residual:rr,alpha:null,beta:null}];
  if(rr<tol)return{x,iterations:0,converged:true,residual:rr,history,message:'El vector inicial ya cumple la tolerancia.'};
  for(let k=1;k<=maxIter;k++){
    const Ap=matVec(A,p),den=dot(p,Ap);
    if(!Number.isFinite(den)||den<=0||!Number.isFinite(rz))return{x,iterations:k-1,converged:false,residual:rr,history,message:'Ruptura numérica en pᵀAp o rᵀz.'};
    const alpha=rz/den;
    if(!Number.isFinite(alpha))return{x,iterations:k-1,converged:false,residual:rr,history,message:'Paso α no finito.'};
    x=x.map((v,i)=>v+alpha*p[i]);
    r=r.map((v,i)=>v-alpha*Ap[i]);
    rr=relativeResidual(A,b,x); // residuo real b-Ax, no solo el recursivo
    const entry:Iteration={k,x:x.slice(),residual:rr,alpha,beta:null};history.push(entry);
    if(!Number.isFinite(rr))return{x,iterations:k,converged:false,residual:rr,history,message:'Apareció un residuo no finito.'};
    if(rr<tol)return{x,iterations:k,converged:true,residual:rr,history,message:'Convergencia alcanzada.'};
    z=apply(r);const next=dot(r,z);
    if(!Number.isFinite(next) || Math.abs(rz)<1e-290)return{x,iterations:k,converged:false,residual:rr,history,message:'Desgaste numérico: producto rᵀz demasiado pequeño.'};
    const beta=next/rz;entry.beta=beta;
    p=z.map((v,i)=>v+beta*p[i]);rz=next;
  }
  return{x,iterations:maxIter,converged:false,residual:rr,history,message:'Se alcanzó el máximo de iteraciones sin cumplir la tolerancia.'};
}
export function parseNumber(value:string):number{
  const s=value.trim().replace(',','.');
  if(!s || !/^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$/.test(s))return NaN;
  return Number(s);
}
export function calcAnalysis(scenario:FinancialScenario,pre:Preconditioner,x0Raw:string[],tolRaw:string,maxRaw:string,omegaRaw:string):Analysis{
  const n=scenario.matrix.length;
  const base:Analysis={valid:false,message:'',symmetric:false,spd:false,n,A:[],b:[],minors:[],eigenvalues:[],eigenPre:[],cond:null,condPre:null,condSym:null,reference:[],preMatrix:[],initialResidual:[],cg:null,pcg:null,errors:[],referenceResidual:null,preconditioner:pre};
  if(n<4||scenario.matrix.some(row=>row.length!==n)||scenario.rhs.length!==n){base.message='El sistema debe ser cuadrado, con n ≥ 4 y un vector b de tamaño n.';return base;}
  const A=scenario.matrix.map(row=>row.map(parseNumber)),b=scenario.rhs.map(parseNumber),x0=x0Raw.map(parseNumber);
  base.A=A;base.b=b;
  if(!isFiniteMatrix(A)||!b.every(Number.isFinite)||x0.length!==n||!x0.every(Number.isFinite)){
    base.message='Hay campos vacíos o no numéricos en A, b o x₀. Completa los valores antes de calcular.';return base;
  }
  const tol=parseNumber(tolRaw),maxIter=parseNumber(maxRaw),omega=parseNumber(omegaRaw);
  if(!Number.isFinite(tol)||tol<=0||tol>=1){base.message='La tolerancia debe estar entre 0 y 1 (valor propuesto: 1e-6).';return base;}
  if(!Number.isInteger(maxIter)||maxIter<1||maxIter>10000){base.message='N máx. debe ser un entero entre 1 y 10.000.';return base;}
  if(pre==='ssor' && (!Number.isFinite(omega)||omega<=0||omega>=2)){base.message='SSOR exige 0 < ω < 2.';return base;}
  base.symmetric=A.every((row,i)=>row.every((v,j)=>
    Math.abs(v-A[j][i])<=1e-10*Math.max(1,Math.abs(v),Math.abs(A[j][i]))));
  if(!base.symmetric){base.message='Matriz no simétrica: A[i,j] debe coincidir con A[j,i]. GCP y GC requieren A simétrica definida positiva.';return base;}
  const fact=cholesky(A);
  if(!fact){base.message='La matriz es simétrica, pero no supera Cholesky: no se demostró que sea definida positiva.';return base;}
  base.spd=true;base.minors=fact.minors;
  try{
    base.reference=solveCholeskyFactor(fact.L,b);
    base.referenceResidual=relativeResidual(A,b,base.reference);
    const raw=conditionSym(A);base.eigenvalues=raw.eig;base.cond=raw.cond;
    const M=precondMatrix(A,pre,pre==='ssor'?omega:1);
    base.preMatrix=M;
    base.initialResidual=sub(b,matVec(A,x0));
    const factorM=cholesky(M);
    if(!factorM)throw Error('No se pudo factorizar el precondicionador M.');
    // Similitud espectral: M⁻¹A ~ L_M⁻¹ A L_M⁻ᵀ, que es simétrica SPD.
    const Linv=inverseLower(factorM.L);
    const B=matMul(matMul(Linv,A),transpose(Linv));
    const nice=B.map((row,i)=>row.map((v,j)=>(v+B[j][i])/2));
    const adjusted=conditionSym(nice);base.eigenPre=adjusted.eig;base.condSym=adjusted.cond;
    // Condición en norma 2 solicitada para M⁻¹A, que no tiene por qué ser simétrica.
    // κ₂(M⁻¹A)=||M⁻¹A||₂ · ||A⁻¹M||₂. Se evalúan las dos normas
    // con los máximos autovalores de las matrices de Gram correspondientes.
    const inverseSPD=(L:Matrix):Matrix => {
      const size=L.length, inv=Array.from({length:size},()=>zeros(size));
      for(let j=0;j<size;j++){
        const e=zeros(size);e[j]=1;
        const col=solveCholeskyFactor(L,e);
        for(let i=0;i<size;i++)inv[i][j]=col[i];
      }
      return inv;
    };
    const MA=matMul(inverseSPD(factorM.L),A), inverseMA=matMul(inverseSPD(fact.L),M);
    const normMA=Math.sqrt(Math.max(0,...eigenSym(matMul(transpose(MA),MA))));
    const normInverseMA=Math.sqrt(Math.max(0,...eigenSym(matMul(transpose(inverseMA),inverseMA))));
    base.condPre=normMA*normInverseMA;
    base.cg=runCG(A,b,x0,tol,maxIter,null);
    base.pcg=runCG(A,b,x0,tol,maxIter,pre==='none'?null:factorM.L);
    base.errors=base.reference.map((value,i):ErrorRow=>{
      const abs=Math.abs(value-base.pcg!.x[i]);
      const relative=Math.abs(value)>1e-14?abs/Math.abs(value):null;
      return {label:`x${i+1}`,reference:value,pcg:base.pcg!.x[i],abs,relative,percent:relative===null?null:relative*100};
    });
    base.valid=true;base.message='Matriz simétrica y definida positiva. Métodos ejecutados correctamente.'+(norm(b)===0?' Atención: b=0, por lo que se usó residuo absoluto para el criterio de parada.':'');
  }catch(e){base.message=e instanceof Error?e.message:'Error numérico durante el análisis.';}
  return base;
}

const UNITS=['Tecnología e infraestructura','Operaciones y personal','Logística y distribución','Comercial y atención'];
export function variableNames(n:number):string[]{return Array.from({length:n},(_,i)=>UNITS[i]??`Proyecto / unidad ${i+1}`);}
function presetNumeric(key:ScenarioKey,n:number):{A:Matrix,b:number[],desired:number[]} {
  const base=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?5.5+((i*3)%4)*0.18:Math.abs(i-j)===1?-.75:0));
  const desired=Array.from({length:n},(_,i)=>[140,120,100,90][i]??Math.max(45,95-i*5));
  if(key==='stress'){
    for(let i=0;i<n;i++)base[i][i]+=0.25;
    desired[0]*=1.12;desired[1]*=1.08;desired[2]*=1.15;desired[3]*=.9;
  }
  let A=base;
  if(key==='ill'){
    const D=Array.from({length:n},(_,i)=>Math.pow(10,(i-(n-1)/2)*3/(n-1)));
    A=base.map((row,i)=>row.map((v,j)=>v*D[i]*D[j]));
  }
  const b=matVec(A,desired);
  return {A,b,desired};
}
export function scenarioPreset(key:ScenarioKey,n:number):FinancialScenario{
  const meta={
    ideal:['Ideal','Planificación estable y coeficientes equilibrados.',
      'Los proyectos presentan penalizaciones similares por desvíos presupuestarios y un acoplamiento moderado entre unidades. A es simétrica, SPD y bien condicionada.',
      'Asignaciones orientativas de referencia en miles de Bs; se mantienen prioridades y costos relativos estables.'],
    stress:['Estrés financiero','Incremento de costos, presión presupuestaria y cambio de prioridades.',
      'Se incrementan penalizaciones marginales en la diagonal y cambia la demanda objetivo: crecen tecnología, operaciones y logística; disminuye comercial. La estructura de interdependencias se conserva.',
      'Refleja un ajuste de asignaciones ante mayores necesidades operativas, inflación y restricciones comerciales.'],
    ill:['Mal condicionado','Desequilibrio severo de escalas entre unidades.',
      'La matriz se forma como D C D, con C SPD y D diagonal positiva de diferentes órdenes de magnitud: conserva simetría/definición positiva, pero aumenta mucho κ(A). Jacobi corrige buena parte del desbalance.',
      'Simula grandes diferencias de sensibilidad o unidades financieras entre departamentos; no representa mediciones de una empresa real.']
  } as const;
  const {A,b}=presetNumeric(key,n), m=meta[key];
  return {key,title:m[0],subtitle:m[1],justification:m[2],financialReading:m[3],matrix:A.map(row=>row.map(v=>String(Number(v.toPrecision(13))))),rhs:b.map(v=>String(Number(v.toPrecision(13))))};
}

@Injectable({providedIn:'root'})
export class NumericalService {
  analyze(s:FinancialScenario,pre:Preconditioner,x0:string[],tol:string,max:string,omega:string):Analysis{
    return calcAnalysis(s,pre,x0,tol,max,omega);
  }
  runTests():QAResult[]{
    const n=4,common={x0:Array(n).fill('0'),tol:'1e-6',max:'100',omega:'1'};
    const analyze=(s:FinancialScenario,p:Preconditioner='jacobi')=>calcAnalysis(s,p,common.x0,common.tol,common.max,common.omega);
    const presets=(['ideal','stress','ill'] as ScenarioKey[]).map(k=>analyze(scenarioPreset(k,n)));
    const t1=presets.every(s=>s.valid&&s.cg?.converged&&s.pcg?.converged);
    const nonsym=scenarioPreset('ideal',n); nonsym.matrix[0][1]='2';const invalidA=analyze(nonsym);
    const notSpd=scenarioPreset('ideal',n);notSpd.matrix[0][0]='-1';const invalidB=analyze(notSpd);
    const base=analyze(scenarioPreset('ideal',n));const edited=scenarioPreset('ideal',n);
    edited.rhs[0]=String(Number(edited.rhs[0])+15);
    const updated=analyze(edited);
    const changed=base.valid&&updated.valid&&Math.abs(base.reference[0]-updated.reference[0])>1e-4;
    const ill=presets[2];const faster=Boolean(ill.valid&&ill.pcg?.converged&&ill.cg?.converged&&ill.pcg.iterations<=ill.cg.iterations);
    const x2=scenarioPreset('ideal',n);x2.matrix[0][0]='';const invalidField=analyze(x2);
    return [
      {code:'QA-01',name:'Carga de los tres escenarios',input:'Ideal, Estrés y Mal condicionado; n=4, x₀=0, tol=10⁻⁶',expected:'SPD y convergencia de GC y GCP en los tres.',obtained:presets.map((a,i)=>`${['Ideal','Estrés','Mal'][i]}: SPD=${a.spd}; GC=${a.cg?.iterations??'—'}; GCP=${a.pcg?.iterations??'—'}`).join(' · '),passed:t1,conclusion:t1?'Conforme':'Revisar iteraciones o escala numérica.'},
      {code:'QA-02',name:'Rechazo de matriz inválida',input:'A₁₂≠A₂₁ y, por separado, A₁₁<0.',expected:'Detener el cálculo y mostrar error explícito.',obtained:`No simétrica: ${!invalidA.valid}; no SPD: ${!invalidB.valid}`,passed:!invalidA.valid&&!invalidB.valid,conclusion:!invalidA.valid&&!invalidB.valid?'Conforme':'Fallo de validación.'},
      {code:'QA-03',name:'Edición dinámica de b',input:'Cambiar b₁ en +15, sin modificar los demás términos.',expected:'Actualizar solución numérica de referencia y GCP.',obtained:`x₁ anterior=${base.reference[0]?.toFixed(6)}; nuevo=${updated.reference[0]?.toFixed(6)}`,passed:changed&&Boolean(updated.pcg?.converged),conclusion:changed?'Conforme':'No se observó recálculo.'},
      {code:'QA-04',name:'Beneficio de Jacobi (mal condicionado)',input:'Escenario Mal condicionado, precondicionador Jacobi y GC simple.',expected:'GCP converge con iteraciones ≤ GC.',obtained:`GC: ${ill.cg?.iterations??'—'} iter.; GCP: ${ill.pcg?.iterations??'—'} iter.`,passed:faster,conclusion:faster?'Conforme (aplica a este caso de prueba).':'No se obtuvo la mejora esperada; revisar datos.'},
      {code:'QA-05',name:'Campo numérico vacío',input:'Vaciar A₁₁.',expected:'Error legible y sin bloqueo de la aplicación.',obtained:invalidField.message,passed:!invalidField.valid&&invalidField.message.includes('vacíos'),conclusion:!invalidField.valid?'Conforme':'Error no detectado.'}
    ];
  }
}
