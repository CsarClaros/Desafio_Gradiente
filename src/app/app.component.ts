import { CommonModule } from '@angular/common';
import { ChangeDetectionStrategy, Component } from '@angular/core';
import { Analysis, FinancialScenario, Preconditioner, QAResult, ScenarioKey, SolverResult } from './models/finance.model';
import { NumericalService, parseNumber, scenarioPreset, variableNames } from './services/numerical.service';

@Component({
  selector:'app-root',standalone:true,imports:[CommonModule],
  templateUrl:'./app.component.html',
  styleUrl:'./app.component.css',
  changeDetection:ChangeDetectionStrategy.Default
})
export class App {
  readonly keys:ScenarioKey[]=['ideal','stress','ill'];
  readonly groups=[
    {id:'caso',label:'01 · Caso de estudio'},
    {id:'formulacion',label:'02 · Formulación'},
    {id:'datos',label:'03 · Entrada de datos'},
    {id:'resultados',label:'04 · Resultados'},
    {id:'comparacion',label:'05 · Comparación'},
    {id:'qa',label:'06 · Pruebas QA'},
    {id:'conclusiones',label:'07 · Conclusiones'}
  ];
  n=4;
  active:ScenarioKey='ideal';
  preconditioner:Preconditioner='jacobi';
  tolerance='1e-6';maxIterations='100';omega='1';
  x0:string[]=Array(4).fill('0');
  scenarios:Record<ScenarioKey,FinancialScenario>={
    ideal:scenarioPreset('ideal',4),stress:scenarioPreset('stress',4),ill:scenarioPreset('ill',4)
  };
  analyses:Record<ScenarioKey,Analysis>={} as Record<ScenarioKey,Analysis>;
  qa:QAResult[]=[];
  selectedModal:''|'model'|'algorithm'|'preconditioner'|'qa'='';
  copied=false;
  readonly sourceUrl='https://cmat.umsa.bo/revmat/rbm4/psunagua_gcp.pdf';
  readonly formatSci=(v:number|null|undefined,dec=4):string=>{
    if(v===null||v===undefined||!Number.isFinite(v))return '—';
    if(v!==0&&(Math.abs(v)>=1e5||Math.abs(v)<1e-4))return v.toExponential(dec);
    return v.toLocaleString('es-BO',{maximumFractionDigits:dec,minimumFractionDigits:dec});
  };
  constructor(private readonly numerical:NumericalService){ this.recalculate(); }
  trackIndex(index:number):number{return index;}
  get current():FinancialScenario{return this.scenarios[this.active];}
  get result():Analysis{return this.analyses[this.active];}
  get names():string[]{return variableNames(this.n);}
  get dimensionIndices():number[]{return Array.from({length:this.n},(_,i)=>i);}
  get totalQA():number{return this.qa.filter(x=>x.passed).length;}
  get anyError():boolean{return this.keys.some(k=>!this.analyses[k]?.valid);}
  get badFinancial():boolean{return !!this.result?.pcg?.x.some(x=>x<0);}
  get scenarioChoices():{key:ScenarioKey;name:string;short:string}[]{return [
    {key:'ideal',name:'Ideal',short:'Estable'},
    {key:'stress',name:'Bajo estrés',short:'Perturbación'},
    {key:'ill',name:'Mal condicionado',short:'Escalas'}
  ];}
  recalculate():void{
    for(const key of this.keys){
      this.analyses[key]=this.numerical.analyze(this.scenarios[key],this.preconditioner,this.x0,this.tolerance,this.maxIterations,this.omega);
    }
  }
  changeScenario(key:ScenarioKey):void{this.active=key;}
  updateMatrix(event:Event,i:number,j:number):void{
    this.current.matrix[i][j]=(event.target as HTMLInputElement).value;
    this.recalculate();
  }
  updateB(event:Event,i:number):void{this.current.rhs[i]=(event.target as HTMLInputElement).value;this.recalculate();}
  updateX0(event:Event,i:number):void{this.x0[i]=(event.target as HTMLInputElement).value;this.recalculate();}
  updateTol(event:Event):void{this.tolerance=(event.target as HTMLInputElement).value;this.recalculate();}
  updateMax(event:Event):void{this.maxIterations=(event.target as HTMLInputElement).value;this.recalculate();}
  updateOmega(event:Event):void{this.omega=(event.target as HTMLInputElement).value;this.recalculate();}
  choosePre(event:Event):void{this.preconditioner=(event.target as HTMLSelectElement).value as Preconditioner;this.recalculate();}
  changeDimension(event:Event):void{
    const next=Number((event.target as HTMLSelectElement).value);
    if(!Number.isInteger(next)||next<4||next>10)return;
    this.n=next;this.x0=Array(next).fill('0');this.resetAll();
  }
  resetAll():void{
    for(const key of this.keys)this.scenarios[key]=scenarioPreset(key,this.n);
    this.qa=[];this.recalculate();
  }
  resetActive():void{this.scenarios[this.active]=scenarioPreset(this.active,this.n);this.recalculate();}
  symmetrize():void{
    const a=this.current.matrix;
    for(let i=0;i<this.n;i++)for(let j=i+1;j<this.n;j++){
      const vi=parseNumber(a[i][j]),vj=parseNumber(a[j][i]);
      if(Number.isFinite(vi)&&Number.isFinite(vj)){
        const avg=String(Number(((vi+vj)/2).toPrecision(12)));
        a[i][j]=avg;a[j][i]=avg;
      }
    }
    this.recalculate();
  }
  runQA():void{this.qa=this.numerical.runTests();}
  calcSpeedup(a:Analysis):string{
    if(!a.cg||!a.pcg)return '—';
    const diff=a.cg.iterations-a.pcg.iterations;
    return diff>0?`${diff} iter. menos`:diff===0?'Mismo número de iter.':`${Math.abs(diff)} iter. más`;
  }
  vectorCompact(v:number[]):string{return v.length?'['+v.map(n=>this.formatSci(n,2)).join(', ')+']':'—';}
  equationText(i:number):string {
    const parts=this.current.matrix[i].map((v,j)=>{
      const x=parseNumber(v);if(!Number.isFinite(x))return `${v}·x${j+1}`;
      const sign=j===0?(x<0?'-':''):(x<0?' − ':' + ');
      return `${sign}${this.formatSci(Math.abs(x),3)}·x${j+1}`;
    });
    return `${parts.join('')} = ${this.current.rhs[i]}`;
  }

  percent(n:number|null|undefined):string{return this.formatSci(n,5)+(n===null||n===undefined?'':' %');}
  chartPoints(s:SolverResult|null, logMin=-12):string{
    if(!s||s.history.length<2)return '';
    const maxK=Math.max(1,this.result.cg?.iterations??0,this.result.pcg?.iterations??0);
    const clamp=(v:number)=>Math.max(logMin,Math.min(0,Math.log10(Math.max(1e-15,v))));
    return s.history.map(v=>`${(46+(v.k/maxK)*666).toFixed(2)},${(44+(clamp(v.residual)/logMin)*190).toFixed(2)}`).join(' ');
  }
  chartEnd(s:SolverResult|null):{x:number,y:number}|null{
    if(!s||!s.history.length)return null;
    const maxK=Math.max(1,this.result.cg?.iterations??0,this.result.pcg?.iterations??0);
    const k=s.history[s.history.length-1];
    return{x:46+k.k/maxK*666,y:44+(Math.max(-12,Math.min(0,Math.log10(Math.max(1e-15,k.residual))))/-12)*190};
  }
  chartTicks():number[]{return [0,-3,-6,-9,-12];}
  chartTolY():number {
    const tol=parseNumber(this.tolerance);
    if(!Number.isFinite(tol)||tol<=0)return 139;
    return 44+Math.max(0,Math.min(1,-Math.log10(tol)/12))*190;
  }
  chartMaxK():number{return Math.max(1,this.result.cg?.iterations??0,this.result.pcg?.iterations??0);}
  errorValues(k:ScenarioKey,type:'abs'|'percent'):number[]{
    return this.analyses[k].errors.map(e=>type==='abs'?e.abs:e.percent??NaN);
  }

  dismissModal():void{this.selectedModal='';}
  averageError(a:Analysis):number|null{return a.valid&&a.errors.length?Math.max(...a.errors.map(e=>e.abs)):null;}
  isNearZero(n:number):boolean{return Math.abs(n)<1e-10;}
  downloadReport():void{
    const lines:string[]=['REPORTE DE RESULTADOS - GCP','COM-254 Métodos Numéricos','Cesar Raymundo Claros Choquemisa',`Fecha de descarga: ${new Date().toLocaleString('es-BO')}`,'',`Precondicionador: ${this.preconditioner}; tolerancia ${this.tolerance}; N máximo ${this.maxIterations}`, ''];
    for(const key of this.keys){
      const a=this.analyses[key]; lines.push(`ESCENARIO: ${this.scenarios[key].title}`,this.scenarios[key].justification,`Estado: ${a.message}`);
      if(a.valid){
        lines.push(`Referencia: ${a.reference.join(', ')}`,`GCP: ${a.pcg?.x.join(', ')}`,`Iteraciones GC: ${a.cg?.iterations}; GCP: ${a.pcg?.iterations}`,`kappa(A): ${a.cond}; kappa(M^-1 A): ${a.condPre}`,`Residuo relativo GCP: ${a.pcg?.residual}`);
      }lines.push('');
    }
    if(this.qa.length){lines.push('PRUEBAS QA');for(const t of this.qa)lines.push(`${t.code} ${t.name}: ${t.passed?'APROBADA':'FALLIDA'} - ${t.obtained}`);}
    const blob=new Blob([lines.join('\n')],{type:'text/plain;charset=utf-8'});
    const url=URL.createObjectURL(blob);const a=document.createElement('a');a.href=url;a.download='reporte-gcp-finanzas.txt';a.click();URL.revokeObjectURL(url);
  }
}
