# Matriz de trazabilidad de la consigna — GCP Finanzas

Base: *Desafío académico de Métodos Numéricos: Aplicación del Método de Gradientes Conjugados Precondicionado en problemas financieros* (cinco páginas, consigna UMSA).

| Código | Exigencia | Implementación | Comprobación |
|---|---|---|---|
| R01 | Individual, caso de finanzas, n≥4 | Caso de planificación para Ingeniería de Sistemas, n editable 4–10 | Secciones Inicio, Caso y Entrada |
| R02 | A simétrica definida positiva | Simetría numérica + Cholesky + menores Sylvester | Módulo resultados y QA-02 |
| R03 | Semántica de variables y coeficientes | x presupuestos en miles Bs, A sensibilidades, b prioridades, contexto | Caso y Formulación |
| R04 | Tres escenarios consistentes | Ideal, estrés y mal escalamiento comparten variables y objetivo | Presets editables y su justificación |
| R05 | Método directo | Cholesky doble sustitución | Solución de referencia en cada escenario |
| R06 | GCP artículo Suñagua | Residuo r0, z=M^-1r, p, alpha, beta y criterio relativo | Servicio numerics + detalle de iteraciones |
| R07 | GC para contraste | Misma A,b,x0,tolerancia sin precondicionar | Tabla de escenario y curvas |
| R08 | Precondicionamiento | Jacobi, SSOR, Identidad; M mostrado numéricamente | Selector + detalle M |
| R09 | Exactitud, errores | Absoluto, relativo (%), relativo global | Tabla por variables y comparativa |
| R10 | κ(A), κ(M^-1A) | κ₂ norma inducida de A y M^-1A, más κ espectral simétrica | Métricas, tabla, README |
| R11 | Reactividad y entrada de datos | Buffer de texto, coma decimal, notación científica, recalculo inmediato | Matriz, b, x0, tol, Nmax, ω |
| R12 | Métodos con invalidaciones | Rechazo de matriz no simétrica, no SPD, entradas inválidas, ω fuera (0,2) | Mensajes y QA-02/QA-05 |
| R13 | Matriz comparativa | Comparación 3 escenarios con soluciones, errores, iter., condición, residuo, conclusión | Tabla obligatoria |
| R14 | QA documentado y ejecutable | Tabla Entrada/Esperado/Obtenido/Conclusión | QA-01 a QA-05 |
| R15 | Bibliografía obligatoria | Suñagua 2020 + Golub/Van Loan + Saad | Conclusiones y README |
| R16 | Repositorio público y web | README, fuente, build y vercel.json | Publicar mediante GitHub/Vercel por estudiante |
| R17 | Sin PDF separado | Toda la teoría integrada al sitio | Cumple Nota de Entregables |
| R18 | Aporte específico Ingeniería de Sistemas | Proceso de decisión y reglas de validación financiera | Caso de estudio |

## Observaciones de interpretación

- La rúbrica hace referencia a “requisitos adicionales estipulados para IA, IS o DS en la Sección 9”, pero en el documento recibido la sección 9 es **Pruebas de software**. Se atendió el énfasis de Ingeniería de Sistemas de la **sección 4**, que está claramente especificado.
- Se usa **n≥4**, y una limitación de interfaz **n≤10** por rendimiento visual/numerical; no impide cumplir el mínimo.
- El documento denomina a veces la solución directa “solución exacta”; siguiendo la propia aclaración terminológica en la sección 6.2, se utiliza «solución numérica de referencia».
- La condición `κ₂(M⁻¹A)` y la condición espectral de la transformación simétrica para PCG **pueden diferir**. Se informan ambas y se indica el significado de cada una.
- La propiedad `GCP≤GC` en iteraciones **no se garantiza universalmente**. QA-04 verifica un caso controlado y reproducible con mal escalamiento diagonal. SSOR se ofrece como alternativa, no se afirma que siempre gane.
- No se afirma que los presets representen datos financieros observados en La Paz. Sus montos son sintéticos, con interpretación económica y modelo documentado.
- No se implementa una asignación estricta de un total presupuestario fijo: Ax=b modela el óptimo de un funcional cuadrático sin restricciones. Una asignación con límites requeriría programación cuadrática restringida.
- Al ser un frontend estático sin persistencia, las ediciones desaparecen si se recarga la página; se incluye exportación de reporte `.txt`, pero no backend ni almacenamiento de datos financieros.
- La entrega real requiere publicar en una URL pública y enlazar un repositorio propio; ninguna URL personal ha sido inventada.

## QA comprobada al producir el proyecto

Se ejecutó un harness de TypeScript aislado del build Angular para evaluar la implementación del servicio. **5/5 pruebas funcionales** aprobadas y **36/36 combinaciones** dimensión/escenario/precondicionador convergieron. Se verificó asimismo el tipado estricto de archivos `.ts` contra definiciones temporales mínimas de Angular; **no equivale** a confirmar la compilación de plantillas por Angular. La instalación npm en este entorno no completó por tiempo de espera de red, por lo que se debe ejecutar `npm install`, `npm run build` y la QA visual en la máquina destino.
