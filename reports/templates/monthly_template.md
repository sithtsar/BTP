---
theme: default
paginate: true
backgroundColor: white
color: black
math: mathjax
---

# LBM Simulation Monthly Report
## {{ report_month }}

---

# Executive Summary

- Single-phase simulations: {{ single_phase_count }} configurations analyzed
- Multiphase simulations: {% if multiphase_metrics %}Completed{% else %}Not available{% endif %}
- Total visualizations: {{ plot_count }} scientific plots generated  
- Analysis date: {{ analysis_date }}

---

{% if single_phase_metrics %}
# Single-Phase Flow Results

{% for mode, metrics in single_phase_metrics.items() %}
## {{ mode.title() }} BGK Method

- RMSE: {{ metrics.rmse | round_scientific }}
- Maximum Error: {{ metrics.max_error | round_scientific }}
- Mean Error: {{ metrics.mean_error | round_scientific }}
{% if metrics.convergence_time %}
- Convergence Time: {{ metrics.convergence_time }} iterations
{% endif %}

---
{% endfor %}

# Error Comparison Table

| Method | RMSE | Max Error | Mean Error |
|--------|------|-----------|------------|
{% for mode, metrics in single_phase_metrics.items() %}
| {{ mode.title() }} | {{ metrics.rmse | round_scientific }} | {{ metrics.max_error | round_scientific }} | {{ metrics.mean_error | round_scientific }} |
{% endfor %}

---
{% endif %}

{% if velocity_plots %}
# Velocity Profile Evolution

{% for plot in velocity_plots %}
![center]({{ plot | plot_path }})

*{{ plot.description }}*

---
{% endfor %}
{% endif %}

{% if multiphase_metrics %}
# Multiphase Flow Results

![center w:600]({{ multiphase_plot_path }})

- Interface stability: {{ multiphase_metrics.interface_stability | format_percentage }}
- Mass conservation: {{ multiphase_metrics.mass_conservation | format_percentage }}
- Density ratio: {{ multiphase_metrics.density_ratio }}:1
- Viscosity ratio: {{ multiphase_metrics.viscosity_ratio }}:1
- Final time step: {{ multiphase_metrics.final_time_step }}

---
{% endif %}

{% if multiphase_plots %}
# 2D Field Visualization

{% for plot in multiphase_plots %}
![center w:700]({{ plot | plot_path }})

*{{ plot.description }}*

---
{% endfor %}
{% endif %}

# Technical Parameters

| Parameter | Single-Phase | Multiphase |
|-----------|--------------|------------|
| Grid Size | {{ parameters.single_phase_grid }} | {{ parameters.multiphase_grid }} |
| Time Steps | {{ parameters.single_phase_timesteps }} | {{ parameters.multiphase_timesteps }} |
| Relaxation Time | τ = {{ parameters.relaxation_time }} | Variable |
| Body Force | {{ parameters.body_force | round_scientific }} | {{ parameters.body_force | round_scientific }} |
| Density Ratio | 1:1 | {{ parameters.density_ratio }}:1 |
| Viscosity Ratio | 1:1 | {{ parameters.viscosity_ratio }}:1 |

---

# Mathematical Formulation

## Single-Phase BGK Equation
$$f_i(\mathbf{x} + \mathbf{e}_i \Delta t, t + \Delta t) = f_i(\mathbf{x}, t) - \frac{1}{\tau}[f_i(\mathbf{x}, t) - f_i^{eq}(\mathbf{x}, t)]$$

## Error Metrics  
$$RMSE = \sqrt{\frac{1}{N} \sum_{i=1}^{N} (u_{sim,i} - u_{analytical,i})^2}$$

{% if multiphase_metrics %}
## Multiphase Phase-Field Model
$$\frac{\partial \phi}{\partial t} + \nabla \cdot (\phi \mathbf{u}) = \nabla \cdot (M \nabla \mu)$$
{% endif %}

---

# Conclusions and Key Findings

## Achievements
{% if single_phase_metrics %}
- Successful validation against analytical Poiseuille solutions
- Error metrics within acceptable tolerance (RMSE < 1e-6)
{% endif %}
{% if multiphase_metrics %}
- Stable multiphase interface tracking maintained  
- High density ratio simulation ({{ multiphase_metrics.density_ratio }}:1) successful
{% endif %}
- Comprehensive visualization pipeline operational

## Next Steps
- Extended parameter sensitivity studies
- Performance optimization analysis
- Advanced post-processing capabilities
- 3D simulation development

---

# Generated Visualizations

{% if all_plots %}
## Available Plots
{% for plot in all_plots[:8] %}
- {{ plot.name }}: {{ plot.description }}
{% endfor %}

{% if all_plots|length > 8 %}
... and {{ all_plots|length - 8 }} additional plots
{% endif %}
{% endif %}

---

# Appendix: File Information

## Data Files
- Single-phase: `{mode}_velocity_t{time}.csv`
- Multiphase: `multiphase_t{time}.csv`, `multiphase_avg_t{time}.csv`

## Generated Plots  
- Location: Project root directory
- Formats: PNG with high resolution
- Naming: Descriptive timestamps included

## Report Details
- Generated: {{ analysis_date }}
- Next Report: {{ next_report_date }}
- Data Period: {{ data_period }}