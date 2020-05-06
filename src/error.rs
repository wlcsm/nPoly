// Taken inspiration from rustyline's error module
use std::fmt;
use std::io;
use std::num;

#[derive(Debug)]
pub enum PolyErr {
    Io(io::Error),
    EmptyPoly,
    Parse(num::ParseIntError),
    ParsePolyError,
}

impl fmt::Display for PolyErr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            PolyErr::Io(ref err) => err.fmt(f),
            PolyErr::EmptyPoly => {
                write!(f, "{}", "Vector to initialise polynomial cannot be empty",)
            }
            PolyErr::Parse(ref err) => err.fmt(f),
            PolyErr::ParsePolyError => write!(f, "{}", "You did it wrong"),
        }
    }
}

impl From<io::Error> for PolyErr {
    fn from(err: io::Error) -> Self {
        PolyErr::Io(err)
    }
}

impl From<num::ParseIntError> for PolyErr {
    fn from(err: num::ParseIntError) -> Self {
        PolyErr::Parse(err)
    }
}

impl From<io::ErrorKind> for PolyErr {
    fn from(kind: io::ErrorKind) -> Self {
        PolyErr::Io(io::Error::from(kind))
    }
}
