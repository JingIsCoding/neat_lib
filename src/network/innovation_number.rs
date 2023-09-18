use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT_INTEGER: AtomicUsize = AtomicUsize::new(1);

pub fn reset() {
    NEXT_INTEGER.store(1, Ordering::SeqCst);
}

pub fn next() -> usize {
    return NEXT_INTEGER.fetch_add(1, Ordering::SeqCst);
}
