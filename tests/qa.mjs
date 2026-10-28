// Pruebas numéricas autocontenidas, sin framework de test adicional.
// Se transpila solamente el servicio TS y se simula @Injectable de Angular.
import { createRequire } from 'node:module';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { resolve, dirname } from 'node:path';
import { strict as assert } from 'node:assert';

const require = createRequire(import.meta.url);
const ts = require('typescript');
const source = resolve(dirname(fileURLToPath(import.meta.url)), '../src/app/services/numerical.service.ts');
const js = ts.transpileModule(readFileSync(source, 'utf8'), {
  compilerOptions: {module:ts.ModuleKind.CommonJS,target:ts.ScriptTarget.ES2022,experimentalDecorators:true},
  reportDiagnostics:true
});
assert.equal(js.diagnostics?.length??0,0,'Errores sintácticos TS');
const exportsMock = {};
const fakeRequire = (name) => {
  if(name==='@angular/core')return {Injectable:()=>(target)=>target};
  throw new Error(`Dependencia inesperada: ${name}`);
};
new Function('require','exports',js.outputText)(fakeRequire,exportsMock);
const {NumericalService,scenarioPreset,calcAnalysis} = exportsMock;
const qa = new NumericalService().runTests();
for(const t of qa){
  console.log(`${t.passed?'APROBADO':'FALLIDO'} ${t.code}: ${t.name}`);
  assert(t.passed,`${t.code}: ${t.obtained}`);
}
for(const n of [4,5,8,10]){
  for(const name of ['ideal','stress','ill']){
    const scenario = scenarioPreset(name,n);
    for(const pre of ['jacobi','ssor','none']){
      const r=calcAnalysis(scenario,pre,Array(n).fill('0'),'1e-6','100','1');
      assert(r.valid&&r.spd,`n=${n}, ${name}, ${pre}: ${r.message}`);
      assert(r.pcg?.converged,`Sin convergencia: n=${n}, ${name}, ${pre}`);
      assert(r.cond>0&&r.condPre>0,'Números de condición no positivos');
      assert.equal(r.errors.length,n,'Errores por variable incompletos');
    }
  }
}
console.log('QA numérica completa: 5 casos funcionales + 36 combinaciones dimensión/escenario/precondicionador.');
