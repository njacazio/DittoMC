# Ditto {#mainpage}

**A fast Monte Carlo generator for minimum-bias events.**

Ditto learns event and particle distributions from a reference Monte Carlo
generator and turns them into a lightweight model designed for extremely fast
event generation.

---

## Start here

### Tune

Build a Ditto model from a reference Monte Carlo sample.

- [Tuning](tuning.md)


### Generate

Produce minimum-bias events from an existing Ditto tune.

- [Generation](generation.md)

### Validate

Compare Ditto with the reference generator and experimental observables.

- [Validation](validation.md)

### API

The API documentation is generated automatically from the C++ source code.

---

## Workflow

```text
Reference Monte Carlo
        │
        ▼
   ┌──────────┐
   │  Tuning  │
   └────┬─────┘
        │
        ▼
   Ditto tune
        │
        ▼
 ┌────────────┐
 │ Generation │
 └─────┬──────┘
       │
       ▼
   Fast events
       │
       ▼
 ┌────────────┐
 │ Validation │
 └────────────┘
```