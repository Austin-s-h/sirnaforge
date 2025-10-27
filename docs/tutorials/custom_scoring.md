# Custom Scoring Tutorial

This tutorial shows how to develop and integrate custom scoring functions for specialized siRNA applications.

## Understanding Thermodynamic Asymmetry Scoring

### What is Thermodynamic Asymmetry?

The thermodynamic asymmetry score is a critical component in siRNA prediction tools that measures the stability difference between the 5' and 3' ends of the guide and passenger strands of the siRNA duplex. This scoring mechanism is based on the fundamental principle that siRNAs function more effectively when the guide strand (antisense) is preferentially incorporated into the RNA-induced silencing complex (RISC).

### The Science Behind Asymmetry Scoring

The asymmetry score predicts which strand of the siRNA duplex will be preferentially loaded into RISC. According to thermodynamic principles:

- The **5' end of the guide strand** should be less thermodynamically stable than the **3' end**
- This asymmetry promotes correct orientation in RISC
- Proper strand selection is crucial for successful RNA interference

The underlying mechanism works because RISC preferentially incorporates the strand whose 5' end is less stable, leading to unwinding of the duplex from the less stable end.

### Research Foundation

The importance of thermodynamic asymmetry in siRNA design is supported by extensive research:

#### Key Publications

1. **Khvorova A et al. (2003)** - "Functional siRNAs and miRNAs exhibit strand bias"
   - Early identification of thermodynamic asymmetry importance for guide strand selection into RISC
   - Demonstrated that RISC incorporation is mediated by thermodynamic asymmetry at position 1 of the duplex
   - Showed that the less stable 5′ nucleotide is preferentially selected

2. **Naito Y et al. (2009)** - "Thermodynamic Stability and Watson-Crick Base Pairing in the Seed Duplex Are Major Determinants of the Efficiency of the siRNA-Based Off-Target Effect"
   - Strong evidence for the relevance of thermodynamic stability in siRNA efficacy
   - Underscores the importance of considering thermodynamic parameters in predicting siRNA efficiency

3. **Amarzguioui M and Prydz H (2004)** - "Designing siRNA That Distinguish between Genes That Differ by a Single Nucleotide"
   - Emphasized the relevance of thermodynamic asymmetry between sense and antisense siRNA strands
   - Identified asymmetry as a critical parameter in designing effective siRNAs

4. **Ichihara M et al. (2017)** - "Simple principles for predicting siRNA efficacy"
   - Comprehensive criteria for effectively predicting siRNA efficacy
   - Includes thermodynamic asymmetry as a key predictive feature

### Implementation in siRNAforge

In siRNAforge, thermodynamic asymmetry scoring is implemented as part of the comprehensive scoring system:

```python
class ThermodynamicAsymmetryScorer:
    """Calculate thermodynamic asymmetry score for siRNA duplex"""

    def calculate_asymmetry_score(self, guide_sequence: str, passenger_sequence: str) -> float:
        """
        Calculate asymmetry score based on 5' vs 3' end stability differences

        Args:
            guide_sequence: Guide strand sequence (antisense)
            passenger_sequence: Passenger strand sequence (sense)

        Returns:
            Asymmetry score (0-1, higher indicates better asymmetry)
        """
        # Calculate 5' end stability (positions 1-4)
        guide_5_prime = guide_sequence[:4]
        passenger_5_prime = passenger_sequence[:4]

        # Calculate 3' end stability (positions -4 to -1)
        guide_3_prime = guide_sequence[-4:]
        passenger_3_prime = passenger_sequence[-4:]

        # Calculate ΔG for each end
        delta_g_5_prime = self._calculate_duplex_energy(guide_5_prime, passenger_5_prime)
        delta_g_3_prime = self._calculate_duplex_energy(guide_3_prime, passenger_3_prime)

        # Asymmetry favors less stable 5' end
        asymmetry = delta_g_3_prime - delta_g_5_prime

        # Normalize to 0-1 scale
        return self._normalize_asymmetry_score(asymmetry)
```

### Scoring Components

The thermodynamic asymmetry score in siRNAforge incorporates several key measurements:

1. **5' End Stability**: Free energy of the first 4 nucleotides of the duplex
2. **3' End Stability**: Free energy of the last 4 nucleotides of the duplex
3. **Asymmetry Ratio**: Difference in stability between 3' and 5' ends
4. **Strand Bias Prediction**: Likelihood of correct guide strand selection

### Comprehensive Thermodynamic Metrics in siRNAforge

siRNAforge evaluates multiple thermodynamic parameters to predict siRNA efficacy. Each metric contributes to the overall assessment of siRNA quality:

#### **GC Content** (`gc_content`)
- **Optimal Range**: 35% to 60% (ideally 40-55%)
- **Rationale**: Balance between duplex stability (enhanced by G:C base pairs) and target accessibility (enabled by A:U base pairs)
- **Too Low (<35%)**: Insufficient duplex stability, poor RISC loading
- **Too High (>60%)**: Overly stable duplex, reduced target accessibility
- **Implementation**: Calculated as percentage of G and C nucleotides in guide sequence

#### **Asymmetry Score** (`asymmetry_score`)
- **Optimal Range**: 0.65 to 1.0 (higher is better)
- **Definition**: Measure of thermodynamic stability difference between 5' and 3' ends of siRNA duplex
- **Mechanism**: Higher asymmetry directs guide strand into RISC and prevents passenger strand loading
- **Calculation**: Based on ΔG difference between duplex ends (ΔG₃' - ΔG₅')
- **Research Basis**: Khvorova et al. (2003) demonstrated this is critical for proper strand selection

#### **Paired Fraction** (`paired_fraction`)
- **Optimal Range**: 0.4 to 0.8 (moderate pairing preferred)
- **Definition**: Fraction of nucleotides in guide strand paired with passenger strand in predicted secondary structure
- **Considerations**:
  - High values (>0.8): Overly stable duplex, reduced strand separation in RISC
  - Low values (<0.4): Insufficient duplex stability for effective loading
  - Moderate values: Optimal balance for RISC processing

#### **Minimum Free Energy (MFE)** (`mfe`)
- **Optimal Range**: -2 to -8 kcal/mol (moderate stability)
- **Definition**: Free energy of the most stable RNA secondary structure
- **Interpretation**:
  - Too negative (<-10 kcal/mol): Overly stable, resistant to RISC processing
  - Too positive (>0 kcal/mol): Unstable duplex, poor RNAi activity
  - Moderate values: Optimal for effective strand separation and target binding

#### **Duplex Stability** (`duplex_stability_dg`)
- **Optimal Range**: -15 to -25 kcal/mol
- **Definition**: Total free energy change (ΔG) for siRNA duplex formation
- **Balance Required**:
  - Too stable (very negative ΔG): Difficult RISC loading and strand separation
  - Too unstable (less negative ΔG): Poor duplex formation and activity
- **Normalized Score**: `duplex_stability_score` converts ΔG to 0-1 scale

#### **End Stability Differences**
- **5' End ΔG** (`dg_5p`): Free energy at guide strand 5' end
- **3' End ΔG** (`dg_3p`): Free energy at guide strand 3' end
- **End Difference** (`delta_dg_end`): Stability difference (ΔG₃' - ΔG₅')

**Optimal Characteristics**:
- **Positive `delta_dg_end`**: 5' end less stable than 3' end (preferred)
- **Target Range**: +2 to +6 kcal/mol difference
- **Mechanism**: Promotes correct guide strand orientation in RISC

#### **Melting Temperature** (`melting_temp_c`)
- **Optimal Range**: 60-78°C for human cells
- **Definition**: Temperature at which 50% of duplexes dissociate
- **Considerations**:
  - High Tm (>80°C): Overly stable duplexes resist strand separation
  - Low Tm (<50°C): Unstable duplexes with poor activity
  - System-specific: Optimal values vary by cell type and organism

#### **Off-target Metrics**
- **Off-target Count** (`off_target_count`): Number of potential off-target sites
- **Transcript Hits** (`transcript_hit_count`, `transcript_hit_fraction`): Target specificity measures
- **Goal**: Minimize off-target binding while maintaining on-target efficacy
- **Acceptable Range**: <5 high-confidence off-targets for most applications

### Interpreting Asymmetry Scores

- **Score 0.8-1.0**: Excellent asymmetry, strong guide strand bias
- **Score 0.6-0.8**: Good asymmetry, likely correct strand selection
- **Score 0.4-0.6**: Moderate asymmetry, may have mixed strand loading
- **Score 0.0-0.4**: Poor asymmetry, risk of passenger strand incorporation

### Integrated Scoring Guidelines

When evaluating siRNA candidates, consider the combined profile of all thermodynamic metrics:

#### **High-Quality siRNA Profile**
```
gc_content: 40-55%
asymmetry_score: 0.7-1.0
paired_fraction: 0.5-0.7
mfe: -4 to -7 kcal/mol
duplex_stability_dg: -18 to -22 kcal/mol
delta_dg_end: +2 to +5 kcal/mol
melting_temp_c: 65-75°C
off_target_count: 0-3
```

#### **Acceptable siRNA Profile**
```
gc_content: 35-60%
asymmetry_score: 0.6-0.7
paired_fraction: 0.4-0.8
mfe: -2 to -8 kcal/mol
duplex_stability_dg: -15 to -25 kcal/mol
delta_dg_end: +1 to +6 kcal/mol
melting_temp_c: 60-78°C
off_target_count: 0-5
```

#### **Poor Quality Indicators**
- GC content <30% or >65%
- Asymmetry score <0.5
- Paired fraction >0.9 (overly stable) or <0.3 (unstable)
- MFE <-10 kcal/mol (too stable) or >0 kcal/mol (unstable)
- Delta_dg_end <0 (wrong asymmetry direction)
- High off-target count (>10)

### Context-Dependent Optimization

**Important Note**: Optimal values can vary based on:

- **Cell Type**: Different cellular environments affect RNA stability
- **Organism**: Species-specific RNA processing differences
- **Delivery Method**: Transfection efficiency affects required stability
- **Target Gene**: Secondary structure and accessibility vary
- **Experimental Conditions**: Temperature, buffer conditions, etc.

#### **Cell Type Considerations**
```python
# Example: Adjust parameters for different cell types
def get_cell_specific_thresholds(cell_type: str) -> dict:
    """Get optimized thresholds for specific cell types"""
    thresholds = {
        'human_hela': {
            'gc_min': 40, 'gc_max': 55,
            'asymmetry_min': 0.7,
            'tm_min': 65, 'tm_max': 75
        },
        'mouse_3t3': {
            'gc_min': 35, 'gc_max': 60,
            'asymmetry_min': 0.65,
            'tm_min': 60, 'tm_max': 78
        },
        'plant_arabidopsis': {
            'gc_min': 45, 'gc_max': 65,
            'asymmetry_min': 0.6,
            'tm_min': 55, 'tm_max': 70
        }
    }
    return thresholds.get(cell_type, thresholds['human_hela'])
```

## Custom Scoring Framework

### Understanding the Scoring System

siRNAforge uses a modular scoring framework where multiple components contribute to the final composite score:

```python
class CompositeScorer:
    """Combines multiple scoring components"""

    def __init__(self):
        self.scorers = {
            'thermodynamic_asymmetry': ThermodynamicAsymmetryScorer(),
            'gc_content': GCContentScorer(),
            'secondary_structure': SecondaryStructureScorer(),
            'off_target': OffTargetScorer()
        }

        self.weights = {
            'thermodynamic_asymmetry': 0.3,  # 30% of total score
            'gc_content': 0.2,               # 20% of total score
            'secondary_structure': 0.2,      # 20% of total score
            'off_target': 0.3                # 30% of total score
        }
```

### Implementing Custom Scorers

You can create custom scoring functions by extending the base scorer class:

```python
from sirnaforge.core.scoring import BaseScorer
from sirnaforge.models.sirna import SiRNACandidate

class CustomThermodynamicScorer(BaseScorer):
    """Custom thermodynamic scoring with specialized parameters"""

    def calculate_score(self, candidate: SiRNACandidate) -> float:
        """
        Calculate custom thermodynamic score

        Args:
            candidate: siRNA candidate to score

        Returns:
            Score between 0-1 (higher is better)
        """
        # Access sequence information
        guide_seq = candidate.guide_sequence
        passenger_seq = candidate.passenger_sequence

        # Custom asymmetry calculation
        asymmetry_score = self._calculate_custom_asymmetry(guide_seq, passenger_seq)

        # Additional thermodynamic factors
        melting_temp_score = self._score_melting_temperature(guide_seq)
        duplex_stability_score = self._score_duplex_stability(guide_seq, passenger_seq)

        # Combine components with custom weights
        final_score = (
            0.5 * asymmetry_score +
            0.3 * melting_temp_score +
            0.2 * duplex_stability_score
        )

        return final_score

    def _calculate_custom_asymmetry(self, guide: str, passenger: str) -> float:
        """Implement custom asymmetry calculation"""
        # Your custom algorithm here
        pass
```

### Integrating Machine Learning Models

Advanced users can integrate machine learning models for sophisticated scoring:

```python
import joblib
from sklearn.ensemble import RandomForestRegressor

class MLThermodynamicScorer(BaseScorer):
    """Machine learning-based thermodynamic scoring"""

    def __init__(self, model_path: str):
        """Load pre-trained ML model"""
        self.model = joblib.load(model_path)

    def calculate_score(self, candidate: SiRNACandidate) -> float:
        """Score using ML model"""
        # Extract features
        features = self._extract_thermodynamic_features(candidate)

        # Predict using ML model
        score = self.model.predict([features])[0]

        # Ensure score is in valid range
        return max(0.0, min(1.0, score))

    def _extract_thermodynamic_features(self, candidate: SiRNACandidate) -> list:
        """Extract features for ML model"""
        features = []

        # Basic sequence features
        features.extend([
            candidate.gc_content,
            len(candidate.guide_sequence),
            candidate.position
        ])

        # Thermodynamic features
        features.extend([
            self._calculate_base_asymmetry(candidate),
            self._calculate_duplex_stability(candidate),
            self._calculate_end_stabilities(candidate)
        ])

        return features
```

## Validation and Performance

### Validating Custom Scores

Always validate your custom scoring functions:

```python
def validate_scorer(scorer, test_candidates):
    """Validate custom scorer performance"""
    scores = []
    for candidate in test_candidates:
        score = scorer.calculate_score(candidate)

        # Check score validity
        assert 0.0 <= score <= 1.0, f"Invalid score: {score}"
        scores.append(score)

    # Check score distribution
    import numpy as np
    print(f"Score range: {np.min(scores):.3f} - {np.max(scores):.3f}")
    print(f"Score mean: {np.mean(scores):.3f}")
    print(f"Score std: {np.std(scores):.3f}")
```

### Performance Optimization

For large-scale analysis, optimize your scoring functions:

```python
from functools import lru_cache

class OptimizedAsymmetryScorer(BaseScorer):
    """Performance-optimized asymmetry scorer"""

    @lru_cache(maxsize=10000)
    def _calculate_duplex_energy(self, seq1: str, seq2: str) -> float:
        """Cached energy calculation"""
        # Expensive calculation cached for repeated sequences
        return self._compute_energy(seq1, seq2)

    def calculate_score_batch(self, candidates: list) -> list:
        """Batch processing for better performance"""
        # Process multiple candidates efficiently
        return [self.calculate_score(c) for c in candidates]
```

## Advanced Topics

### Research Integration

When developing scoring functions, always reference current research:

1. **Literature Review**: Stay current with siRNA efficacy research
2. **Experimental Validation**: Test against known effective siRNAs
3. **Benchmarking**: Compare against established scoring methods
4. **Publication**: Consider publishing novel scoring approaches

### Future Directions

The field of siRNA scoring continues to evolve:

- **Deep Learning Models**: Neural networks for complex pattern recognition
- **Multi-Modal Scoring**: Integration of sequence, structure, and expression data
- **Species-Specific Models**: Tailored scoring for different organisms
- **Real-Time Optimization**: Adaptive scoring based on experimental feedback

For more advanced implementations, refer to:
- {doc}`../api_reference` for complete API documentation
- {doc}`../developer/architecture` for scoring system design
- `src/sirnaforge/core/design.py` for current scoring implementation
