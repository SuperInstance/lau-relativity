//! Agent spacetime modeling — relativistic analogs for agent reference frames.

use serde::{Serialize, Deserialize};
use crate::lorentz::gamma;
use crate::minkowski::FourVector;
use crate::kinematics::ReferenceFrame;

/// An agent's "spacetime position" in a multi-agent system.
/// Uses relativistic analogies: agents at different "velocities" (processing rates)
/// experience subjective time differently.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgentSpacetime {
    /// Agent identifier.
    pub agent_id: String,
    /// Subjective proper time (time experienced by the agent).
    pub proper_time: f64,
    /// Coordinate time (wall clock / system time).
    pub coordinate_time: f64,
    /// "Velocity" — ratio of the agent's processing speed to the reference rate.
    /// beta = 0 means at rest (same as reference), beta → 1 means very fast.
    pub beta: f64,
    /// Position in some abstract state space.
    pub state_position: [f64; 3],
}

impl AgentSpacetime {
    pub fn new(agent_id: &str, coordinate_time: f64, beta: f64, state_position: [f64; 3]) -> Self {
        let proper_time = coordinate_time / gamma(beta);
        Self {
            agent_id: agent_id.to_string(),
            proper_time,
            coordinate_time,
            beta,
            state_position,
        }
    }

    /// How much subjective time passes for each unit of coordinate time.
    pub fn time_dilation_factor(&self) -> f64 {
        1.0 / gamma(self.beta)
    }

    /// Advance the agent by dt coordinate time.
    pub fn advance(&mut self, dt: f64) {
        self.coordinate_time += dt;
        self.proper_time += dt / gamma(self.beta);
    }

    /// Lorentz contract a distance as seen from this agent's frame.
    pub fn contracted_distance(&self, distance: f64) -> f64 {
        distance / gamma(self.beta)
    }

    /// The "spacetime interval" between this agent and another.
    /// s² = (Δt)² - (Δx)²/c² where Δx is state-space distance and c is a normalization.
    pub fn interval_to(&self, other: &AgentSpacetime) -> f64 {
        let dt = other.coordinate_time - self.coordinate_time;
        let dx = self.state_distance(other);
        let c_norm = 1.0; // Normalized c=1 in state space
        dt * dt - dx * dx / (c_norm * c_norm)
    }

    /// Euclidean distance in state space.
    pub fn state_distance(&self, other: &AgentSpacetime) -> f64 {
        let d0 = self.state_position[0] - other.state_position[0];
        let d1 = self.state_position[1] - other.state_position[1];
        let d2 = self.state_position[2] - other.state_position[2];
        (d0 * d0 + d1 * d1 + d2 * d2).sqrt()
    }

    /// Whether this agent can causally influence the other.
    pub fn can_influence(&self, other: &AgentSpacetime) -> bool {
        self.interval_to(other) >= 0.0 // timelike or lightlike
    }
}

/// A multi-agent reference frame system.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgentFrameSystem {
    pub agents: Vec<AgentSpacetime>,
    /// Reference frame velocity.
    pub reference_beta: f64,
}

impl AgentFrameSystem {
    pub fn new(reference_beta: f64) -> Self {
        Self { agents: Vec::new(), reference_beta }
    }

    /// Add an agent to the system.
    pub fn add_agent(&mut self, agent: AgentSpacetime) {
        self.agents.push(agent);
    }

    /// Transform all agents to a new reference frame with given beta.
    /// Agents' coordinate times and positions are Lorentz-transformed.
    pub fn boost_to_frame(&self, new_beta: f64) -> Vec<AgentSpacetime> {
        let relative_beta = (new_beta - self.reference_beta) / (1.0 - new_beta * self.reference_beta);
        let g = gamma(relative_beta);
        self.agents.iter().map(|a| {
            let new_coord_time = g * (a.coordinate_time - relative_beta * a.state_position[0]);
            let new_x = g * (a.state_position[0] - relative_beta * a.coordinate_time);
            AgentSpacetime {
                agent_id: a.agent_id.clone(),
                proper_time: a.proper_time, // proper time is invariant
                coordinate_time: new_coord_time,
                beta: relative_beta,
                state_position: [new_x, a.state_position[1], a.state_position[2]],
            }
        }).collect()
    }

    /// Find all agents that can causally influence a given agent.
    pub fn causal_past(&self, target_idx: usize) -> Vec<usize> {
        if target_idx >= self.agents.len() {
            return Vec::new();
        }
        let target = &self.agents[target_idx];
        self.agents.iter().enumerate()
            .filter(|(i, a)| *i != target_idx && a.can_influence(target))
            .map(|(i, _)| i)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_agent_proper_time_dilation() {
        let agent = AgentSpacetime::new("a1", 10.0, 0.6, [0.0, 0.0, 0.0]);
        let g = gamma(0.6);
        assert!((agent.proper_time - 10.0 / g).abs() < 1e-10);
    }

    #[test]
    fn test_agent_advance() {
        let mut agent = AgentSpacetime::new("a1", 0.0, 0.0, [0.0, 0.0, 0.0]);
        agent.advance(5.0);
        assert!((agent.coordinate_time - 5.0).abs() < 1e-10);
        assert!((agent.proper_time - 5.0).abs() < 1e-10); // at rest
    }

    #[test]
    fn test_agent_advance_moving() {
        let mut agent = AgentSpacetime::new("a1", 0.0, 0.8, [0.0, 0.0, 0.0]);
        agent.advance(10.0);
        assert!((agent.coordinate_time - 10.0).abs() < 1e-10);
        let g = gamma(0.8);
        assert!((agent.proper_time - 10.0 / g).abs() < 1e-8);
    }

    #[test]
    fn test_causal_influence_timelike() {
        let a = AgentSpacetime::new("a1", 0.0, 0.0, [0.0, 0.0, 0.0]);
        let b = AgentSpacetime::new("a2", 10.0, 0.0, [5.0, 0.0, 0.0]);
        assert!(a.can_influence(&b)); // timelike separation (100 - 25 > 0)
    }

    #[test]
    fn test_no_causal_influence_spacelike() {
        let a = AgentSpacetime::new("a1", 0.0, 0.0, [0.0, 0.0, 0.0]);
        let b = AgentSpacetime::new("a2", 1.0, 0.0, [100.0, 0.0, 0.0]);
        assert!(!a.can_influence(&b)); // spacelike separation
    }

    #[test]
    fn test_frame_system_boost() {
        let mut system = AgentFrameSystem::new(0.0);
        system.add_agent(AgentSpacetime::new("a1", 10.0, 0.0, [5.0, 0.0, 0.0]));
        let boosted = system.boost_to_frame(0.0);
        assert_eq!(boosted.len(), 1);
        // Boost with beta=0 should be identity
        assert!((boosted[0].coordinate_time - 10.0).abs() < 1e-8);
    }

    #[test]
    fn test_contracted_distance() {
        let agent = AgentSpacetime::new("a1", 0.0, 0.6, [0.0, 0.0, 0.0]);
        let d = agent.contracted_distance(10.0);
        let g = gamma(0.6);
        assert!((d - 10.0 / g).abs() < 1e-10);
    }

    #[test]
    fn test_causal_past() {
        let mut system = AgentFrameSystem::new(0.0);
        system.add_agent(AgentSpacetime::new("source", 0.0, 0.0, [0.0, 0.0, 0.0]));
        system.add_agent(AgentSpacetime::new("reachable", 5.0, 0.0, [1.0, 0.0, 0.0]));
        system.add_agent(AgentSpacetime::new("unreachable", 1.0, 0.0, [100.0, 0.0, 0.0]));
        let past = system.causal_past(1);
        assert!(past.contains(&0));
        assert!(!past.contains(&2));
    }

    #[test]
    fn test_state_distance() {
        let a = AgentSpacetime::new("a1", 0.0, 0.0, [1.0, 2.0, 2.0]);
        let b = AgentSpacetime::new("a2", 0.0, 0.0, [0.0, 0.0, 0.0]);
        assert!((a.state_distance(&b) - 3.0).abs() < 1e-10);
    }
}
