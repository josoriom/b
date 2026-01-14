pub mod decode;
pub use decode::decode;
pub mod decode2;
pub use decode2::decode2;
pub mod encode2;
pub use encode2::encode2;
pub mod encode;
pub use encode::encode;
pub mod utilities;

#[cfg(test)]
mod tests;
