# siRNAforge Architecture

Short version (TL;DR)
---------------------

- What it is: A layered Python toolkit for end-to-end siRNA design.
- How it’s organized: CLI → Workflow → Core → Models → Data → Pipeline, with shared Utils/Validation.
- How data flows: Gene query → transcripts → design → scores → ranked output.
- How to extend: plug in new scorers, data providers, and output writers via small interfaces.

Use this page as a map. Skim the diagrams and the “Architectural Layers” section first; dive into details only as needed.

## Overview

siRNAforge is built as a modern Python package with clear separation of concerns, type safety, and extensibility. The architecture follows domain-driven design principles with distinct layers for different responsibilities.

```{raw} html
<div class="mermaid">
graph TD
    A[CLI Interface] --> B[Workflow Orchestration]
    B --> C[Core Algorithms]
    B --> D[Data Access Layer]
    B --> E[Pipeline Integration]

    C --> C1[siRNA Design]
    C --> C2[Thermodynamics]
    C --> C3[Off-target Analysis]

    D --> D1[Gene Search]
    D --> D2[ORF Analysis]
    D --> D3[External APIs]

    E --> E1[Nextflow Pipeline]
    E --> E2[Docker Integration]

    F[Data Models] --> C
    F --> D
    G[Utils & Validation] --> C
    G --> D

    style A fill:#e1f5fe
    style B fill:#f3e5f5
    style C fill:#e8f5e8
    style D fill:#fff3e0
    style E fill:#fce4ec
</div>
```

## System Architecture Flow

```{raw} html
<div class="mermaid">
sequenceDiagram
    participant User
    participant CLI
    participant Workflow
    participant Core
    participant Data
    participant External

    User->>CLI: sirnaforge workflow GENE
    CLI->>Workflow: orchestrate_workflow(gene)
    Workflow->>Data: search_gene_transcripts(gene)
    Data->>External: query_ncbi/ensembl()
    External-->>Data: transcript_sequences
    Data-->>Workflow: validated_transcripts

    Workflow->>Core: design_sirnas(transcripts)
    Core->>Core: thermodynamic_analysis()
    Core->>Core: off_target_prediction()
    Core-->>Workflow: scored_candidates

    Workflow->>CLI: results_summary
    CLI->>User: formatted_output + files
</div>
```

## Package Structure

```{raw} html
<div class="mermaid">
graph TB
    subgraph "siRNAforge Package"
        A[__init__.py<br/>Version & Config]
        B[cli.py<br/>Command Interface]
        C[workflow.py<br/>Orchestration]

        subgraph "Core Algorithms"
            D[design.py<br/>siRNA Design]
            E[thermodynamics.py<br/>RNA Folding]
            F[off_target.py<br/>Off-target Analysis]
        end

        subgraph "Data Models"
            G[sirna.py<br/>Pydantic Models]
            H[schemas.py<br/>Pandera Schemas]
        end

        subgraph "Data Access"
            I[gene_search.py<br/>Gene/Transcript Search]
            J[orf_analysis.py<br/>ORF Analysis]
            K[base.py<br/>Base Classes]
        end

        subgraph "Pipeline Integration"
            L[nextflow/<br/>Pipeline Modules]
            M[resources/<br/>Resource Mgmt]
        end

        subgraph "Utilities & Validation"
            N[utils/<br/>Logging Utils]
            O[validation/<br/>Validation System]
        end
    end

    B --> C
    C --> D
    C --> E
    C --> F
    C --> I
    C --> J

    D --> G
    E --> G
    F --> G
    I --> G
    J --> G

    C --> L

    style D fill:#e8f5e8
    style E fill:#e8f5e8
    style F fill:#e8f5e8
    style G fill:#fff3e0
    style H fill:#fff3e0
    style I fill:#f3e5f5
    style J fill:#f3e5f5
    style K fill:#f3e5f5
</div>
```

### Directory Structure
```
src/sirnaforge/
├── __init__.py          # Package initialization and version
├── cli.py              # Command-line interface (Typer/Rich)
├── workflow.py         # High-level workflow orchestration
│
├── core/               # Core algorithms and business logic
│   ├── design.py      # siRNA design algorithms
│   ├── thermodynamics.py  # RNA folding and energy calculations
│   └── off_target.py  # Off-target prediction algorithms
│
├── models/             # Data models and validation
│   ├── sirna.py       # Pydantic models for siRNA data
│   └── schemas.py     # Pandera validation schemas
│
├── data/               # Data access and external APIs
│   ├── base.py        # Base classes for data providers
│   ├── gene_search.py # Gene/transcript search functionality
│   └── orf_analysis.py # Open reading frame analysis
│
├── pipeline/           # Pipeline and workflow integration
│   ├── nextflow/      # Nextflow workflow configs and runners
│   └── resources/     # Resource management
│
├── utils/              # Shared utilities
│   └── logging_utils.py # Logging configuration
│
└── validation/         # Data validation and QC
    ├── config.py      # Validation configuration
    ├── middleware.py  # Validation middleware
    └── utils.py       # Validation utilities
```

## Architectural Layers

### 1. CLI Layer (`cli.py`)

**Purpose**: User interface and command orchestration

**Technologies**:
- [Typer](https://typer.tiangolo.com/) for CLI framework
- [Rich](https://rich.readthedocs.io/) for beautiful output

**Responsibilities**:
- Command parsing and validation
- Progress indicators and user feedback
- Error handling and user-friendly messages
- Configuration management

### Configuration & Default Resolution (`sirnaforge.config`)

The `sirnaforge.config.reference_policy` module centralizes how user inputs and default references are resolved. `WorkflowInputSpec` captures raw CLI/API inputs, while `ReferencePolicyResolver` produces metadata-rich `ReferenceChoice` objects that indicate whether a reference was explicitly provided, auto-selected, or disabled. The workflow records these choices (currently for transcriptome references) inside `logs/workflow_summary.json` so production runs can be audited without inspecting CLI arguments.

### 2. Workflow Layer (`workflow.py`)

**Purpose**: High-level process orchestration

**Pattern**: Facade/Coordinator

**Responsibilities**:
- Coordinate multi-step workflows
- Handle data flow between components
- Manage temporary files and outputs
- Provide consistent APIs for different entry points

### 3. Core Layer (`core/`)

**Purpose**: Core algorithms and business logic

**Pattern**: Strategy/Template Method

#### Core Components:

##### `design.py` - siRNA Design Engine
```python
class SiRNADesigner:
    """Main siRNA design orchestrator"""

class SiRNACandidate:
    """Individual siRNA candidate with scoring"""

class DesignParameters:
    """Configuration for design algorithms"""
```

##### `thermodynamics.py` - RNA Structure Analysis
```python
class ThermodynamicsCalculator:
    """ViennaRNA integration for structure prediction and asymmetry scoring"""

class ThermodynamicAsymmetryScorer:
    """Calculate thermodynamic asymmetry for guide strand selection"""

class SecondaryStructure:
    """RNA secondary structure representation"""
```

**Thermodynamic Asymmetry Implementation:**

The thermodynamic asymmetry scoring is a critical component that predicts guide strand selection into RISC (RNA-induced silencing complex). This implementation is based on research showing that siRNAs with less stable 5' ends on the guide strand are more effectively incorporated into RISC.

**Key Research Foundation:**
- Khvorova A et al. (2003): Demonstrated thermodynamic asymmetry importance for RISC incorporation
- Naito Y et al. (2009): Established thermodynamic stability as a major determinant of siRNA efficiency
- Amarzguioui M and Prydz H (2004): Identified asymmetry as critical for distinguishing target genes
- Ichihara M et al. (2017): Comprehensive principles including thermodynamic asymmetry for efficacy prediction

**Algorithm Components:**
1. **5' End Stability Analysis**: Calculates free energy of duplex 5' terminus (positions 1-4)
2. **3' End Stability Analysis**: Calculates free energy of duplex 3' terminus (positions -4 to -1)
3. **Asymmetry Ratio Calculation**: Measures stability difference (ΔG₃' - ΔG₅')
4. **Strand Bias Prediction**: Predicts likelihood of correct guide strand selection

##### `off_target.py` - Specificity Analysis
```python
class OffTargetPredictor:
    """Multi-genome off-target prediction"""

class AlignmentResult:
    """Off-target alignment with scoring"""
```

### 4. Model Layer (`models/`)

**Purpose**: Data models and validation

**Technologies**: [Pydantic](https://pydantic.dev/) for data validation

**Pattern**: Data Transfer Objects (DTOs)

#### Key Models:

```python
class SiRNACandidate(BaseModel):
    """Complete siRNA candidate with all metadata"""
    sirna_id: str
    guide_sequence: str
    passenger_sequence: str
    position: int
    composite_score: float
    asymmetry_score: float
    gc_content: float

class DesignParameters(BaseModel):
    """Design configuration with validation"""
    sirna_length: int = Field(ge=19, le=23)
    gc_min: float = Field(ge=0, le=100)
    gc_max: float = Field(ge=0, le=100)
    top_n: int = Field(ge=1)

class FilterCriteria(BaseModel):
    """Quality filters for candidate selection"""
    gc_min: float = 30.0
    gc_max: float = 60.0
    max_poly_runs: int = 3
    min_asymmetry_score: float = 0.65

class ScoringWeights(BaseModel):
    """Relative weights for composite scoring"""
    asymmetry: float = 0.25
    gc_content: float = 0.20
    accessibility: float = 0.25
    off_target: float = 0.20
    empirical: float = 0.10
```

### 5. Data Layer (`data/`)

**Purpose**: External data access and integration

**Pattern**: Repository/Adapter

#### Data Providers:

##### `gene_search.py` - Gene Information Retrieval
```python
class GeneSearcher:
    """Multi-database gene search"""

class GeneSearchResult:
    """Complete gene search result"""

class EnsemblClient(AbstractDatabaseClient):
    """Ensembl REST API integration"""

class RefSeqClient(AbstractDatabaseClient):
    """RefSeq/NCBI integration"""

class GencodeClient(AbstractDatabaseClient):
    """GENCODE database integration"""
```

##### `orf_analysis.py` - Sequence Analysis
```python
class ORFAnalyzer:
    """Open reading frame validation"""
```

### 6. Validation Layer (`validation/`)

**Purpose**: Data validation and quality control

**Pattern**: Middleware/Decorator

**Technologies**: [Pandera](https://pandera.readthedocs.io/) for schema validation

#### Validation Components:

##### `config.py` - Validation Configuration
```python
class ValidationConfig(BaseModel):
    """Configuration for validation system"""

class ValidationLevel(str, Enum):
    """Validation strictness levels"""
    STRICT = "strict"
    WARNING = "warning"
    DISABLED = "disabled"

class ValidationStage(str, Enum):
    """Pipeline stages for validation"""
    INPUT = "input"
    DESIGN = "design"
    OUTPUT = "output"
```

##### `utils.py` - Validation Utilities
```python
class ValidationResult:
    """Container for validation results"""

def validate_sirna_candidates(df: pd.DataFrame) -> ValidationResult:
    """Validate siRNA candidate data"""

def validate_fasta_sequences(sequences: list) -> ValidationResult:
    """Validate FASTA sequence data"""
```

### 7. Pipeline Layer (`pipeline/`)

**Purpose**: External pipeline integration

**Technologies**: Nextflow, Docker

**Responsibilities**:
- Nextflow workflow orchestration
- Docker container management
- Batch processing coordination
- Resource management

#### Pipeline Components:

##### `nextflow/` - Workflow Management
```python
class NextflowConfig:
    """Nextflow execution configuration"""

class NextflowRunner:
    """Nextflow workflow execution"""
```

##### `resources/` - Resource Management
```python
class ResourceManager:
    """Compute resource allocation and monitoring"""
```

## Design Principles

### 1. Type Safety

All components use comprehensive type hints and Pydantic validation:

```python
from typing import List, Optional
from pydantic import BaseModel, Field

class DesignParameters(BaseModel):
    sirna_length: int = Field(21, ge=19, le=23)
    top_candidates: int = Field(10, ge=1)
    gc_content_range: tuple[float, float] = (30.0, 60.0)
```

### 2. Separation of Concerns

Each layer has distinct responsibilities:
- **CLI**: User interaction
- **Workflow**: Process orchestration
- **Core**: Algorithm implementation
- **Models**: Data representation
- **Data**: External integration

### 3. Dependency Injection

Components use constructor injection and composition patterns:

```python
class SiRNADesigner:
    def __init__(self, parameters: DesignParameters) -> None:
        """Initialize designer with configuration parameters."""
        self.parameters = parameters
        # Components are instantiated as needed
        # ThermodynamicCalculator() instantiated when required

class GeneSearcher:
    def __init__(self, timeout: int = 30, max_retries: int = 3):
        """Initialize with configurable database clients."""
        self.clients: dict[DatabaseType, AbstractDatabaseClient] = {
            DatabaseType.ENSEMBL: EnsemblClient(timeout=timeout),
            DatabaseType.REFSEQ: RefSeqClient(timeout=timeout),
            DatabaseType.GENCODE: GencodeClient(timeout=timeout),
        }
```

### 4. Error Handling

Comprehensive error handling with custom exceptions:

```python
class SiRNAForgeException(Exception):
    """Base exception for all siRNAforge errors"""

class DesignException(SiRNAForgeException):
    """siRNA design specific errors"""

class ValidationException(SiRNAForgeException):
    """Input validation errors"""
```

### 5. Configuration Management

Centralized configuration with environment support:

```python
class WorkflowConfig(BaseModel):
    """Workflow-specific configuration"""

class ValidationConfig(BaseModel):
    """Validation configuration with environment variable support"""

class NextflowConfig(BaseModel):
    """Nextflow execution parameters"""

    class Config:
        env_prefix = "SIRNAFORGE_"
```

## Data Flow

### 1. Complete Workflow

```{mermaid}
graph TD
    A[Gene Query] --> B[Gene Search]
    B --> C[Transcript Retrieval]
    C --> D[ORF Analysis]
    D --> E[siRNA Design]
    E --> F[Thermodynamic Analysis]
    F --> G[Off-target Prediction]
    G --> H[Scoring & Ranking]
    H --> I[Output Generation]
```

### 2. Component Interaction

```{mermaid}
graph LR
    CLI --> Workflow
    Workflow --> Core
    Core --> Models
    Core --> Data
    Pipeline --> Core
    Utils --> All[All Components]
```

## Extension Points

### 1. Custom Scoring Functions

The current implementation uses fixed scoring within `SiRNADesigner._score_candidates()`. Custom scoring can be implemented by:

```python
class CustomSiRNADesigner(SiRNADesigner):
    def _score_candidates(self, candidates: list[SiRNACandidate]) -> list[SiRNACandidate]:
        """Override with custom scoring logic"""
        for candidate in candidates:
            candidate.composite_score = self._custom_score(candidate)
        return candidates

    def _custom_score(self, candidate: SiRNACandidate) -> float:
        # Custom scoring logic
        return score
```

### 2. Additional Data Sources

Extend the database client system:

```python
class CustomDatabaseClient(AbstractDatabaseClient):
    async def search_gene(self, query: str) -> GeneSearchResult:
        # Custom gene search logic
        return result

# Register with the searcher
searcher = GeneSearcher()
searcher.clients[DatabaseType.CUSTOM] = CustomDatabaseClient()
```

### 3. New Validation Rules

Extend the validation system:

```python
class CustomValidationRules:
    def validate_custom_criteria(self, candidates: pd.DataFrame) -> ValidationResult:
        # Custom validation logic
        return result

# Use with validation middleware
from sirnaforge.validation.middleware import ValidationMiddleware
validator = ValidationMiddleware(custom_rules=CustomValidationRules())
```

## Performance Considerations

### 1. Asynchronous Operations

External API calls and I/O operations use asyncio:

```python
async def search_multiple_databases(query: str) -> List[SearchResult]:
    tasks = [
        search_ensembl(query),
        search_refseq(query),
        search_gencode(query)
    ]
    return await asyncio.gather(*tasks)
```

### 2. Memory Management

Large datasets are processed in chunks:

```python
def design_from_large_file(file_path: Path) -> Iterator[SiRNACandidate]:
    for chunk in read_fasta_chunks(file_path, chunk_size=1000):
        yield from design_candidates(chunk)
```

### 3. Caching

Expensive computations are cached:

```python
from functools import lru_cache

@lru_cache(maxsize=1000)
def calculate_thermodynamics(sequence: str) -> ThermodynamicResult:
    # Expensive ViennaRNA calculation
    return result
```

## Testing Architecture

### 1. Unit Tests (`tests/unit/`)
- Test individual components in isolation
- Mock external dependencies
- Focus on algorithm correctness

### 2. Integration Tests (`tests/integration/`)
- Test component interactions
- Use real external services (with rate limiting)
- Validate end-to-end workflows

### 3. Pipeline Tests (`tests/pipeline/`)
- Test Nextflow pipeline components
- Container-based testing
- Resource usage validation

## Deployment Architecture

### 1. Local Development
- `uv` for dependency management
- Direct Python execution
- Local debugging and testing

### 2. Container Deployment
- Multi-stage Docker builds
- Optimized for size and security
- Environment-specific configurations

### 3. Pipeline Deployment
- Nextflow for workflow orchestration
- Support for multiple execution platforms
- Resource management and monitoring

## Future Architecture Considerations

### 1. Microservices
- Potential split into specialized services
- API gateway for service coordination
- Independent scaling of components

### 2. Cloud Integration
- Cloud storage for large datasets
- Serverless functions for lightweight operations
- Managed services for databases

### 3. Plugin System
- Dynamic loading of algorithms
- Third-party extensions
- Community contributions

This architecture provides a solid foundation for current needs while maintaining flexibility for future enhancements.
