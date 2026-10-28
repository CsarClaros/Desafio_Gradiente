# Método de Gradientes Conjugados Precondicionado — Planificación financiera

**Universidad Mayor de San Andrés · COM-254 Métodos Numéricos · Desafío individual**  
**Estudiante:** Cesar Raymundo Claros Choquemisa  
**Docente:** Lic. Brigida Carvajal  
**Mención de enfoque propuesta:** Ingeniería de Sistemas.

Aplicación web académica desarrollada en **Angular 21 + TypeScript 5.9 + Tailwind CSS 3**, basada en el desafío *Aplicación del Método de Gradientes Conjugados Precondicionado en problemas financieros*. Funciona 100 % en el navegador, sin backend, registros o API. Incluye los tres escenarios, solver directo de referencia, GC, GCP, números de condición, comparaciones, visualización de convergencia, documentación y QA en la web.

## 1. Instalación y ejecución

Requisitos: Node.js `>=20.19` (línea 20), `>=22.12` (línea 22) o línea 24 admitida por Angular 21.2; npm funcional.

```bash
npm install
npm start
```

Abre `http://localhost:4200`. Para generar producción:

```bash
npm run build
```

### Paquetes coherentes (sin mezcla Tailwind v3/v4)

| Tecnología | Versión |
|---|---|
| Angular (core, compiler, CLI, build) | `21.2.4` |
| TypeScript | `5.9.3` |
| Tailwind CSS | `3.4.19` |
| PostCSS | `8.5.6` |
| Autoprefixer | `10.4.21` |
| Zone.js | `0.16.3` |

Tailwind 3 funciona con `@tailwind base; @tailwind components; @tailwind utilities;` en `src/styles.css`, `postcss.config.js` y `tailwind.config.js`. Se mantienen las clases académicas `.card`, `.field`, `.data-input`, `.btn-primary`, etc.

## 2. Problema financiero y sentido de los datos

Se simula la asignación de montos expresados en **miles de bolivianos** entre cuatro unidades empresariales (y unidades adicionales si aumentas n):

1. Tecnología e infraestructura.
2. Operaciones y personal.
3. Logística y distribución.
4. Comercial y atención.

El sistema deriva de la minimización de la función convexa

\[
J(x)=\tfrac12 x^\top A x-b^\top x,\quad\nabla J(x)=Ax-b=0.
\]

- **x:** propuesta de asignaciones de recursos, miles de Bs.
- **A (simétrica definida positiva):** matriz de sensibilidades marginales y acoplamiento financiero entre unidades. Las diagonales penalizan los desvíos y los coeficientes fuera de diagonal modelan las interdependencias. Los coeficientes se utilizan en un modelo normalizado.
- **b:** vector de prioridades marginales ponderadas. Los casos de prueba predeterminados se **calibran didácticamente** a partir de asignaciones objetivo `t` creando `b=A·t`, para ofrecer resultados verificables. **No son registros financieros reales** de empresas bolivianas.

Como es un problema cuadrático **sin restricciones**, no impone `sum(x)=B`, `x≥0` ni cotas de inversión. Para distribuir un presupuesto máximo fijo habría que formular un problema de optimización cuadrática sujeto a restricciones y utilizar algoritmos acordes. Aquí el foco matemático exigido por el desafío es resolver correctamente `Ax=b` mediante GCP.

### Tres escenarios

- **Ideal:** matriz tridiagonal SPD con diagonales entre 5.5 y 6.04 y enlaces de `−0.75` entre unidades contiguas. Buena condición espectral. Demandas calibradas de referencia en miles de Bs.
- **Bajo estrés:** aumenta la penalización diagonal en 0.25 y cambian las asignaciones objetivo por presiones en tecnología, operaciones y logística, con ajuste comercial. Se conserva el esquema SPD.
- **Mal condicionado:** se construye `A = D C D`, donde C es la matriz SPD de referencia y D es diagonal positiva con escalas geométricas entre `10⁻¹·⁵` y `10¹·⁵`. Se obtienen sensibilidades muy desbalanceadas, número de condición alrededor de `10^6` para 4 variables y mejora medible al usar Jacobi; la SPD se conserva matemáticamente.

## 3. Requerimientos matemáticos implementados

### Matriz simétrica y definida positiva

Primero se comprueba `|A_ij−A_ji|≤1e-10·escala`. Después se intenta la factorización `A=L Lᵀ` de Cholesky. Si falla, la aplicación detiene los métodos. Se muestran los menores principales `Δ_k=∏_{i≤k} L_ii²`, estrictamente positivos para los casos SPD. Esto implementa el criterio de Sylvester.

### Método directo de referencia

Se resuelven sucesivamente `Ly=b` y `Lᵀx=y`. Esta es una **solución numérica de referencia** en coma flotante; no se presenta como solución simbólica exacta. Se informa el residuo obtenido.

### Gradientes conjugados (GC) y precondicionado (GCP)

- Residuo inicial `r₀=b−Ax₀`.
- Resolver `Mzₖ=rₖ` (o `zₖ=rₖ` para GC).
- `p₀=z₀`.
- `αₖ=(rₖᵀzₖ)/(pₖᵀApₖ)`.
- `xₖ₊₁=xₖ+αₖpₖ`.
- `rₖ₊₁=rₖ−αₖApₖ`.
- `βₖ=(rₖ₊₁ᵀzₖ₊₁)/(rₖᵀzₖ)`.
- `pₖ₊₁=zₖ₊₁+βₖpₖ`.

Se evalúa el **residuo verdadero** `||b−Ax||₂/||b||₂ < 1e−6` por defecto (si `b=0`, se usa residuo absoluto y se informa). Se admiten `x₀`, tolerancia y máximo de iteraciones editables. Se registran las iteraciones, α, β, x y norma residual, y se dibujan las curvas de GC/GCP.

### Precondicionadores

- **Jacobi:** `M=diag(A)`, adecuado para mal escalamiento diagonal.
- **SSOR:** `M=(D+ωL)D⁻¹(D+ωU)/[ω(2−ω)]`, con `0<ω<2` y `A=D+L+U`. Alternativa para acoplamiento; no garantiza mejora general.
- **Sin precondicionar:** `M=I`, lo que reproduce GC.

Se obtiene `κ₂(A)` mediante autovalores de A con rotaciones de Jacobi simétrico. Para el precondicionado se muestran **dos indicadores diferentes**: (1) `κ₂(M⁻¹A) = ||M⁻¹A||₂ · ||A⁻¹M||₂`, calculado usando las mayores singularidades obtenidas a través de matrices de Gram e inversas; (2) la razón de autovalores positivos de la matriz simétrica `L_M⁻¹ A L_M⁻ᵀ` (`M=L_M L_Mᵀ`), que es semejante a `M⁻¹A` y explica la convergencia teórica del GCP. En general **no son iguales**; distinguirlos es esencial.

### Errores

- Absoluto componente: `Eaᵢ=|x_ref,ᵢ−x_GCP,ᵢ|`.
- Relativo: `Erᵢ=Eaᵢ/|x_ref,ᵢ|` si `x_ref,ᵢ≠0`.
- Porcentual: `Epᵢ=100 Erᵢ`.
- Residuo relativo: `||b−Ax_GCP||₂/||b||₂` (caso `b≠0`).

## 4. Arquitectura Angular

```text
src/
├── main.ts
├── index.html
├── styles.css                         # Tailwind 3 y componentes visuales
└── app/
    ├── app.component.ts               # Estado reactivo por escenario, inputs, exportación
    ├── app.component.html             # Módulos académicos y modal
    ├── app.component.css
    ├── models/
    │   └── finance.model.ts           # Tipos fuertes
    └── services/
        └── numerical.service.ts       # Cholesky, PCG, SSOR, autovalores, QA
```

Los **inputs numéricos son de tipo `text` con `inputmode="decimal"`**, con un buffer de texto persistente por escenario. Esto permite introducir valores parciales (`-`, `1e`, `0.`) sin perder foco ni destruir/recrear el input. Las filas utilizan `trackBy` por índice fijo y nunca se reconstruyen durante la edición. La aplicación valida antes de calcular; admite coma decimal y notación científica, pero el usuario no necesita reiniciar el navegador si deja una celda temporalmente incompleta. Los cambios se recalculan al instante y se preservan por escenario.

## 5. Módulo QA

En la propia web, pulsar **«Ejecutar pruebas automáticas»**. Las pruebas emplean **copias** de presets de 4 variables y no alteran el estado editado:

| Prueba | Entrada | Esperado |
|---|---|---|
| QA-01 | Cargar tres escenarios | Simetría, SPD y convergencia de GC y GCP |
| QA-02 | Romper simetría / SPD | Mensaje de error sin cerrar la app |
| QA-03 | Cambiar `b₁` en `+15` | Recalcular referencia y solución GCP |
| QA-04 | Mal condicionado con Jacobi | GCP converge en iteraciones `≤` GC |
| QA-05 | Campo A₁₁ vacío | Validación legible, sin ruptura |

La sección QA documenta **Entrada, Resultado esperado, Resultado obtenido y Conclusión**, sin proclamar aprobadas pruebas que el usuario todavía no haya ejecutado. Hay una prueba **manual** adicional para verificar retención del foco y del texto parcial.

También puede ejecutarse una verificación numérica de consola después de `npm install`:

```bash
npm run test:qa
```

Este comando no sustituye las pruebas de navegación visual en el navegador.

## 6. Resultados verificables de referencia (n = 4; Jacobi; tol = 10⁻⁶)

Aproximados, obtenidos con pruebas aisladas de los algoritmos (no atribuidos a un build Angular compilado):

| Escenario | κ₂(A) | κ₂(M⁻¹A) | κ espectral simétrica | GC iter. | GCP iter. |
|---|---:|---:|---:|---:|
| Ideal | 1.534 | 1.524 | 1.524 | 4 | 4 |
| Estrés | 1.507 | 1.497 | 1.497 | 4 | 4 |
| Mal condicionado | ~1,051,167 | 9.510 | 1.524 | 5 | 4 |

**Precaución:** un método puede alcanzar un residuo pequeño y, pese a ello, tener errores relativos relevantes si el sistema tiene alta condición. Precondicionar busca mejorar la convergencia espectral, no garantizar en todas las matrices menor cantidad de iteraciones.

## 7. Publicación y entregables

La consigna exige **únicamente**: página web pública completamente funcional, sin credenciales, y repositorio público con README e historial de commits. **No requiere informe PDF separado**; toda fundamentación, resultados, guía de uso y QA está integrada en la web.

1. Crear repositorio público y subir `src/`, archivos de configuración, `.gitignore`, README y `package-lock.json` después de la instalación.
2. Ejecutar `npm run build` y corregir cualquier problema de compilación.
3. Publicar en Vercel conectando GitHub (Output directory `dist/gradientes-conjugados-finanzas/browser`) o mediante la carpeta compilada según la plataforma. Verificar públicamente la URL.
4. Ejecutar las 5 pruebas QA y la prueba manual de inputs **en la versión publicada**. Registrar URL de sitio y repositorio para entrega.

## 8. Fuentes

- **Fuente académica exigida:** Suñagua S., P. (2020). *Método de Gradientes Conjugados Precondicionado*. Revista Boliviana de Matemática – UMSA 04, 2–7. [PDF oficial](https://cmat.umsa.bo/revmat/rbm4/psunagua_gcp.pdf).
- Golub, G. H. y Van Loan, C. F. (2013). *Matrix Computations*, 4th ed. Johns Hopkins University Press.
- Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems*, 2nd ed. SIAM.
- [Versiones de Angular (oficial)](https://angular.dev/reference/versions).

**Limitaciones declaradas:** datos financieros sintéticos, optimización sin restricciones, sin respaldo persistente de cambios al recargar la página, y límite de tamaño de interfaz de 4 a 10 variables para conservar legibilidad. Los valores de condición para matrices muy mal condicionadas pueden degradarse por aritmética de coma flotante. La QA comparativa es específica del escenario de prueba, no una ley universal sobre GCP.
